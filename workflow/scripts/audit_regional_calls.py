"""Trace P11 regional calls; this is not a rerun of the full-genome caller.

Use the portable bundle and the review environment (pysam 0.23.3).
Keep windows unchanged: heuristic BAQ can depend on calling context.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pysam
import pysam.bcftools as bcftools

TARGETS = [('I', 27486, 'G', 'A'), ('IV', 1364942, 'C', 'T'),
           ('IX', 37309, 'A', 'G')]
FILTER = 'QUAL="." || QUAL<30 || FMT/DP="." || FMT/DP<8 || FMT/DP>200'


def fingerprint(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_records(path, targets):
    positions = {(chrom, pos) for chrom, pos, _, _ in targets}
    records = []
    with open(path) as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            f = line.rstrip().split('\t')
            if (f[0], int(f[1])) not in positions:
                continue
            if len(f) != 10:
                raise ValueError('Expected one sample')
            info = dict(item.split('=', 1) if '=' in item else (item, True)
                        for item in f[7].split(';'))
            records.append({'chrom': f[0], 'pos_1based': int(f[1]),
                            'ref': f[3], 'alt': f[4],
                            'qual': None if f[5] == '.' else float(f[5]),
                            'filter': f[6], 'info': info,
                            'sample': dict(zip(f[8].split(':'), f[9].split(':'))),
                            'vcf_line': line.rstrip()})
    return records


def rejection_reasons(record):
    reasons = []
    qual = record['qual']
    dp = record['sample'].get('DP', '.')
    if qual is None or qual < 30:
        reasons.append('missing_QUAL' if qual is None else 'QUAL<30')
    if dp == '.':
        reasons.append('missing_FMT/DP')
    elif int(dp) < 8 or int(dp) > 200:
        reasons.append('FMT/DP<8' if int(dp) < 8 else 'FMT/DP>200')
    return reasons


def audit(bundle, output, targets=TARGETS):
    bundle, output = Path(bundle), Path(output)
    if output.exists():
        raise FileExistsError(output)
    required = ['reference.fa', 'reference.fa.fai', 'extraction_windows.bed']
    for tech in ('illumina', 'ont'):
        required += [f'P11.{tech}.regional.bam', f'P11.{tech}.regional.bam.bai']
    for name in required:
        if not (bundle / name).is_file():
            raise FileNotFoundError(bundle / name)
    with pysam.FastaFile(str(bundle / 'reference.fa')) as fasta:
        for chrom, pos, ref, alt in targets:
            if fasta.fetch(chrom, pos - 1, pos) != ref or alt not in 'ACGT':
                raise ValueError('Target/reference mismatch')
    output.mkdir(parents=True)
    result = {
        'scope': 'Regional reconstruction, not original full-genome intermediates',
        'pysam_version': pysam.__version__,
        'samtools_version': pysam.__samtools_version__,
        'script_sha256': fingerprint(__file__),
        'sources_sha256': {name: fingerprint(bundle / name) for name in required},
        'filter_expression': FILTER,
        'targets': [dict(zip(('chrom', 'pos_1based', 'ref', 'alt'), t)) for t in targets],
        'experiments': {},
    }
    experiments = [('illumina_default', 'illumina', 13, []),
                   ('ont_bq7', 'ont', 7, ['-B', '--skip-indels']),
                   ('illumina_no_baq', 'illumina', 13, ['-B'])]
    for name, tech, bq, extra in experiments:
        paths = {s: str(output / f'{name}.{s}.vcf')
                 for s in ('pileup', 'raw', 'normalized', 'pass')}
        calls = [
            ('mpileup', extra + ['-q', '20', '-Q', str(bq), '-d', '1000',
                '--ns', '3844', '-a', 'FORMAT/DP,FORMAT/AD', '-Ov',
                '-o', paths['pileup'], '-f', str(bundle / 'reference.fa'),
                '-R', str(bundle / 'extraction_windows.bed'),
                str(bundle / f'P11.{tech}.regional.bam')]),
            ('call', ['--ploidy', '2', '-mv', '-Ov', '-o', paths['raw'], paths['pileup']]),
            ('norm', ['-f', str(bundle / 'reference.fa'), '-c', 'e', '-m', '-any',
                      '-Ov', '-o', paths['normalized'], paths['raw']]),
            ('view', ['-e', FILTER, '-Ov', '-o', paths['pass'], paths['normalized']]),
        ]
        logs = []
        for command, args in calls:
            dispatch = getattr(bcftools, command)
            dispatch(*args, catch_stdout=False)
            logs.append({'command': command, 'args': args, 'messages': dispatch.get_messages()})
        records = {stage: read_records(path, targets) for stage, path in paths.items()}
        for record in records['normalized']:
            record['rejection_reasons'] = rejection_reasons(record)
            key = (record['chrom'], record['pos_1based'], record['ref'], record['alt'])
            passed = any((r['chrom'], r['pos_1based'], r['ref'], r['alt']) == key
                         for r in records['pass'])
            if passed != (not record['rejection_reasons']):
                raise ValueError('Predicate explanation disagrees with bcftools view')
        header = Path(paths['pileup']).read_text().splitlines()
        versions = [line for line in header if line.startswith('##bcftoolsVersion=')]
        result['experiments'][name] = {'commands': logs, 'bcftools_version_header': versions,
                                       'records': records}
    (output / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    audit(args.bundle, args.output)
