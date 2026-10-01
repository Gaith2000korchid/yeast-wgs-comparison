"""Verify P11 archive reads, select paired Illumina reads, and combine ONT runs.

Pair selection uses SHA256(seed + ':' + normalized read name), independent of
Python's random implementation. Each mate is selected by the same decision.
"""
import argparse
import gzip
import hashlib
import json
import re
from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path

from download_reads import verify


def records(stream):
    while True:
        header = stream.readline()
        if not header:
            return
        seq, plus, qual = stream.readline(), stream.readline(), stream.readline()
        if not header.startswith(b'@') or not plus.startswith(b'+') or not qual:
            raise ValueError('Malformed or truncated four-line FASTQ')
        sequence, quality = seq.rstrip(b'\r\n'), qual.rstrip(b'\r\n')
        if not sequence or len(sequence) != len(quality):
            raise ValueError('FASTQ sequence/quality length mismatch')
        if min(quality) < 33 or max(quality) > 126:
            raise ValueError('Expected Phred+33 FASTQ quality characters')
        yield header, seq, plus, qual


def read_name(record):
    name = record[0].split()[0][1:]
    if name.endswith((b'/1', b'/2')):
        name = name[:-2]
    return name


def select_pair(name, fraction, seed):
    threshold = int(fraction * (1 << 64))
    number = int.from_bytes(hashlib.sha256(str(seed).encode() + b':' + name).digest()[:8], 'big')
    return number < threshold


def subsample_pairs(r1, r2, out1, out2, fraction, seed):
    if not 0 < fraction <= 1:
        raise ValueError('Fraction must be in (0, 1]')
    counts = {'input_pairs': 0, 'input_bases': 0, 'selected_pairs': 0, 'selected_bases': 0,
              'fraction_requested': fraction, 'seed': seed, 'algorithm': 'SHA256 seed:name, first 64 bits'}
    with ExitStack() as stack:
        inputs = [stack.enter_context(gzip.open(p, 'rb')) for p in (r1, r2)]
        outputs = []
        for path in (out1, out2):
            raw = stack.enter_context(open(path, 'wb'))
            outputs.append(stack.enter_context(gzip.GzipFile(filename='', fileobj=raw, mode='wb', compresslevel=1, mtime=0)))
        for left, right in zip_longest(records(inputs[0]), records(inputs[1])):
            if left is None or right is None or read_name(left) != read_name(right):
                raise ValueError('Illumina mates are not synchronized')
            bases = len(left[1].rstrip()) + len(right[1].rstrip())
            counts['input_pairs'] += 1
            counts['input_bases'] += bases
            if select_pair(read_name(left), fraction, seed):
                counts['selected_pairs'] += 1
                counts['selected_bases'] += bases
                for record, stream in zip((left, right), outputs):
                    stream.write(b''.join(record))
    return counts


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def ont_identity(record):
    tokens = record[0].decode().strip().split()
    if len(tokens) < 2:
        raise ValueError('Missing original historical ONT read name')
    name = tokens[1].removesuffix('/1')
    match = re.search(r'_(twodirections|template|complement)_(pass|fail)_', name)
    if not match:
        raise ValueError(f'Unknown historical ONT header: {name}')
    molecule = name[:match.start()] + '_' + name[match.end():]
    return match.group(1), match.group(2), molecule


def combine_ont(paths, output):
    # Keep consensus 2D reads (pass and fail); template/complement are correlated.
    seen = set()
    summaries = []
    with open(output, 'wb') as raw, gzip.GzipFile(filename='', fileobj=raw, mode='wb', compresslevel=1, mtime=0) as dest:
        for path in paths:
            counts = {'filename': Path(path).name, 'archive_reads': 0, 'archive_bases': 0,
                      'selected_2d_reads': 0, 'selected_2d_bases': 0, 'read_types': {}, 'header_examples': []}
            with gzip.open(path, 'rb') as stream:
                for record in records(stream):
                    kind, status, molecule = ont_identity(record)
                    bases = len(record[1].rstrip())
                    counts['archive_reads'] += 1
                    counts['archive_bases'] += bases
                    key = kind + '_' + status
                    counts['read_types'][key] = counts['read_types'].get(key, 0) + 1
                    if kind != 'twodirections':
                        continue
                    if molecule in seen:
                        raise ValueError(f'Duplicate ONT consensus molecule: {molecule}')
                    seen.add(molecule)
                    counts['selected_2d_reads'] += 1
                    counts['selected_2d_bases'] += bases
                    if len(counts['header_examples']) < 3:
                        counts['header_examples'].append(record[0].decode().strip())
                    dest.write(b''.join(record))
            summaries.append(counts)
    if not seen:
        raise ValueError('No ONT consensus 2D reads found')
    return summaries


def prepare(manifest, raw, out, fraction, seed):
    metadata = json.loads(Path(manifest).read_text())
    entries = metadata['files']
    for entry in entries:
        if entry['sample_accession'] != metadata['sample_accession']:
            raise ValueError('Manifest includes an unmatched sample')
        if Path(entry['filename']).name != entry['filename'] or not verify(Path(raw) / entry['filename'], entry):
            raise ValueError(f"Archive validation failed: {entry['filename']}")
    illumina = [e for e in entries if e['platform'] == 'ILLUMINA']
    ont = [e for e in entries if e['platform'] == 'OXFORD_NANOPORE']
    if len(illumina) != 2 or not illumina[0]['filename'].endswith('_1.fastq.gz') or not illumina[1]['filename'].endswith('_2.fastq.gz') or not ont:
        raise ValueError('Expected one ordered Illumina pair and at least one ONT run')
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    paths = [out / 'P11.R1.fastq.gz', out / 'P11.R2.fastq.gz', out / 'P11.ont.fastq.gz']
    partials = [p.with_name(p.name + '.part') for p in paths]
    counts = subsample_pairs(*(Path(raw) / e['filename'] for e in illumina), *partials[:2], fraction, seed)
    ont_stats = combine_ont([Path(raw) / e['filename'] for e in ont], partials[2])
    for partial, target in zip(partials, paths):
        if target.is_file() and sha256(partial) == sha256(target):
            partial.unlink()
        else:
            partial.replace(target)
    provenance = {'sample': metadata['sample'], 'sample_accession': metadata['sample_accession'],
                  'archive_files': entries, 'illumina': counts, 'ont_selection': 'twodirections only, pass and fail; no template/complement', 'ont_runs': ont_stats,
                  'prepared_sha256': {p.name: sha256(p) for p in paths}}
    (out / 'preparation.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps({'illumina': counts, 'ont_reads': sum(s['selected_2d_reads'] for s in ont_stats), 'ont_bases': sum(s['selected_2d_bases'] for s in ont_stats)}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', default='config/P11.downloads.json')
    p.add_argument('--raw', default='data/real/raw')
    p.add_argument('--out', default='data/real/prepared')
    p.add_argument('--fraction', type=float, default=0.25)
    p.add_argument('--seed', type=int, default=20261001)
    a = p.parse_args()
    prepare(a.manifest, a.raw, a.out, a.fraction, a.seed)
