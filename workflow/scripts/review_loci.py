"""Measure diagnostic fragment/read support at an explicit list of SNP loci.

Counts are independent of variant calling: no BAQ, no genotype likelihoods.
Overlapping Illumina mates contribute once; conflicting eligible mates are
excluded rather than used to choose an allele. This is not a truth classifier.
"""
import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import hashlib
import json
import os
from pathlib import Path
from statistics import median
import xml.etree.ElementTree as ET

EXCLUDE_FLAGS = 3844


def read_loci(path):
    with open(path, newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        required = {'chrom', 'pos_1based', 'ref', 'alt', 'category', 'illumina_gt', 'ont_gt'}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError('Locus TSV is missing required columns')
        rows = []
        seen = set()
        for row in reader:
            row['pos_1based'] = int(row['pos_1based'])
            if row['pos_1based'] < 1 or row['ref'] not in 'ACGT' or len(row['ref']) != 1 or row['alt'] not in 'ACGT' or len(row['alt']) != 1 or row['ref'] == row['alt']:
                raise ValueError('Expected positive 1-based coordinates and distinct A/C/G/T SNP alleles')
            if row['category'] not in {'shared', 'illumina_only', 'ont_only'}:
                raise ValueError('Unknown comparison category')
            key = (row['chrom'], row['pos_1based'], row['ref'], row['alt'])
            if key in seen:
                raise ValueError('Duplicate locus')
            seen.add(key)
            rows.append({key: row[key] for key in ['chrom', 'pos_1based', 'ref', 'alt', 'category', 'illumina_gt', 'ont_gt']})
    if not rows:
        raise ValueError('Empty locus list')
    return rows


def observation(read, pos0):
    for query_pos, reference_pos in read.get_aligned_pairs():
        if reference_pos != pos0:
            continue
        if query_pos is None:
            return None, 'deletion_or_refskip'
        if read.query_sequence is None or read.query_qualities is None:
            return None, 'missing_sequence_or_quality'
        base = read.query_sequence[query_pos].upper()
        if base not in 'ACGT':
            base = 'N'
        distance = min(query_pos - read.query_alignment_start, read.query_alignment_end - 1 - query_pos)
        return {'name': read.query_name, 'base': base, 'quality': int(read.query_qualities[query_pos]),
                'reverse': bool(read.is_reverse), 'read1': bool(read.is_read1),
                'aligned_end_distance': distance,
                'soft_clipped': any(op == 4 for op, length in (read.cigartuples or []) if length)}, None
    return None, 'no_base_at_position'


def collect(bam, chrom, pos0, minimum_mapping_quality):
    observations = []
    reasons = Counter()
    for read in bam.fetch(chrom, pos0, pos0 + 1):
        reasons['alignments_overlapping_position'] += 1
        if read.flag & EXCLUDE_FLAGS:
            reasons['excluded_alignment_flags'] += 1
            continue
        if read.mapping_quality < minimum_mapping_quality:
            reasons['below_mapping_quality'] += 1
            continue
        item, reason = observation(read, pos0)
        if reason:
            reasons[reason] += 1
        else:
            observations.append(item)
    return observations, dict(reasons)


def support(observations, ref, alt, minimum_base_quality, paired):
    groups = defaultdict(list)
    for i, item in enumerate(observations):
        if item['quality'] >= minimum_base_quality:
            # Illumina query name identifies the fragment; ONT counts records.
            groups[item['name'] if paired else i].append(item)
    retained = []
    conflicts = collapsed = 0
    for items in groups.values():
        if len({item['base'] for item in items}) > 1:
            conflicts += 1
            continue
        collapsed += len(items) - 1
        # Prefer highest original quality, then read1 on an agreeing-mate tie.
        retained.append(max(items, key=lambda item: (item['quality'], item['read1'], item['aligned_end_distance'])))
    counts = Counter(item['base'] for item in retained)
    alternate = [item for item in retained if item['base'] == alt]
    depth = len(retained)
    return {'min_base_quality': minimum_base_quality, 'depth': depth,
            'ref_count': counts[ref], 'alt_count': counts[alt],
            'other_count': depth - counts[ref] - counts[alt],
            'base_counts': {base: counts[base] for base in 'ACGTN'},
            'alt_fraction': counts[alt] / depth if depth else None,
            'ref_fraction': counts[ref] / depth if depth else None,
            'alt_forward': sum(not item['reverse'] for item in alternate),
            'alt_reverse': sum(item['reverse'] for item in alternate),
            'alt_median_original_base_quality': median(item['quality'] for item in alternate) if alternate else None,
            'alt_within_10bp_of_aligned_query_end': sum(item['aligned_end_distance'] <= 10 for item in alternate),
            'alt_on_soft_clipped_alignment': sum(item['soft_clipped'] for item in alternate),
            'overlapping_mate_observations_collapsed': collapsed,
            'conflicting_eligible_fragments_excluded': conflicts}


def digest(path):
    value = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def igv_session(path, reference, bams, locus):
    parent = Path(path).resolve().parent
    relative = lambda resource: os.path.relpath(Path(resource).resolve(), parent).replace(os.sep, '/')
    root = ET.Element('Session', genome=relative(reference), locus=locus, version='3')
    resources = ET.SubElement(root, 'Resources')
    for label, bam in bams.items():
        ET.SubElement(resources, 'Resource', name=label, path=relative(bam), index=relative(str(bam) + '.bai'))
    ET.indent(root)
    ET.ElementTree(root).write(path, encoding='utf-8', xml_declaration=True)


def review(loci, reference, illumina, ont, output, mapping_quality=20, illumina_bq=13, ont_bqs=(13, 7)):
    import pysam
    loci, reference, illumina, ont = map(str, (loci, reference, illumina, ont))
    if mapping_quality < 0 or illumina_bq < 0 or not ont_bqs or any(q < 0 for q in ont_bqs) or len(set(ont_bqs)) != len(ont_bqs):
        raise ValueError('Quality cutoffs must be distinct nonnegative integers')
    rows = read_loci(loci)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    # Verify every reference allele before writing any evidence/session files.
    with ExitStack() as stack:
        fasta = stack.enter_context(pysam.FastaFile(reference))
        bams = {p: stack.enter_context(pysam.AlignmentFile(path, 'rb')) for p, path in [('illumina', illumina), ('ont', ont)]}
        for row in rows:
            chrom, pos0 = row['chrom'], row['pos_1based'] - 1
            if chrom not in fasta.references or pos0 >= fasta.get_reference_length(chrom):
                raise ValueError(f'Out-of-reference coordinate: {chrom}:{pos0+1}')
            if fasta.fetch(chrom, pos0, pos0 + 1).upper() != row['ref']:
                raise ValueError(f'Reference allele mismatch: {chrom}:{pos0+1}')
            for bam in bams.values():
                if chrom not in bam.references or bam.get_reference_length(chrom) != fasta.get_reference_length(chrom):
                    raise ValueError('BAM/reference contig names or lengths differ')
                if not bam.has_index():
                    raise ValueError('BAM must be indexed')
        evidence = []
        for row in rows:
            result = dict(row)
            result['views'] = {}
            result['alignment_filters'] = {}
            for platform, qualities in [('illumina', [illumina_bq]), ('ont', ont_bqs)]:
                observations, reasons = collect(bams[platform], row['chrom'], row['pos_1based']-1, mapping_quality)
                result['alignment_filters'][platform] = reasons
                for bq in qualities:
                    result['views'][f'{platform}_bq{bq}'] = support(observations, row['ref'], row['alt'], bq, platform == 'illumina')
            evidence.append(result)
        first = rows[0]
        length = fasta.get_reference_length(first['chrom'])
        location = f"{first['chrom']}:{max(1,first['pos_1based']-100)}-{min(length,first['pos_1based']+100)}"
    sources = {'loci': loci, 'reference': reference, 'illumina_bam': illumina, 'ont_bam': ont,
               'reference_index': reference + '.fai', 'illumina_bam_index': illumina + '.bai', 'ont_bam_index': ont + '.bai'}
    metadata = {'pysam_version': pysam.__version__, 'pysam_samtools_version': pysam.__samtools_version__,
                'script_sha256': digest(__file__), 'min_mapping_quality': mapping_quality, 'excluded_flags': EXCLUDE_FLAGS,
                'sources': {key: {'path': str(path), 'sha256': digest(path)} for key, path in sources.items()},
                'interpretation': 'Diagnostic original-base support, not caller DP/BAQ/GQ, genotype likelihoods or truth. Eligible overlapping mates agree or are excluded; agreeing mates count once. Loci are a selected convenience sample, not genome-wide statistics.',
                'visual_review_status': 'not_performed', 'loci': evidence}
    (out/'read_support.json').write_text(json.dumps(metadata, indent=2) + '\n')
    fields = ['chrom', 'pos_1based', 'ref', 'alt', 'category', 'illumina_gt', 'ont_gt', 'view', 'min_base_quality',
              'depth', 'ref_count', 'alt_count', 'other_count', 'alt_fraction', 'alt_forward', 'alt_reverse',
              'alt_median_original_base_quality', 'alt_within_10bp_of_aligned_query_end', 'alt_on_soft_clipped_alignment',
              'overlapping_mate_observations_collapsed', 'conflicting_eligible_fragments_excluded']
    with (out/'read_support.tsv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for item in evidence:
            for view, values in item['views'].items():
                record = {key: item[key] for key in fields[:7]}
                record.update(view=view)
                record.update({key: values[key] for key in fields[8:]})
                writer.writerow(record)
    with (out/'regions.bed').open('w') as stream:
        for row in rows:
            stream.write(f"{row['chrom']}\t{row['pos_1based']-1}\t{row['pos_1based']}\t{row['category']}:{row['ref']}>{row['alt']}\n")
    igv_session(out/'igv_session.xml', reference, {'Illumina': illumina, 'ONT historical consensus 2D': ont}, location)
    with (out/'visual_review_template.tsv').open('w', newline='') as stream:
        fields = ['chrom', 'pos_1based', 'ref', 'alt', 'category', 'reviewer', 'review_date', 'igv_version', 'visual_status', 'notes', 'snapshot_filename']
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({**{key:row[key] for key in fields[:5]}, 'visual_status':'not_reviewed'})
    print(f"Measured {len(rows)} selected loci; visual review remains pending. Outputs: {out}")


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--loci', default='docs/results/P11/bq7.loci_for_review.tsv')
    p.add_argument('--reference', default='data/real/reference.fa')
    p.add_argument('--illumina', default='results/P11/bam/P11.illumina.bam')
    p.add_argument('--ont', default='results/P11/bam/P11.ont.bam')
    p.add_argument('--out', default='results/P11/review')
    p.add_argument('--min-mapping-quality', type=int, default=20)
    p.add_argument('--illumina-base-quality', type=int, default=13)
    p.add_argument('--ont-base-qualities', nargs='+', type=int, default=[13, 7])
    a = p.parse_args()
    review(a.loci, a.reference, a.illumina, a.ont, a.out, a.min_mapping_quality, a.illumina_base_quality, tuple(a.ont_base_qualities))
