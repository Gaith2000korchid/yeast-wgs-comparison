"""Exact normalized alternate-SNP overlap in a shared depth mask; not accuracy."""
import argparse
import bisect
import csv
import gzip
import json
from collections import defaultdict


def read_bed(path):
    intervals = defaultdict(list)
    with open(path) as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip().split("\t")
            start, end = int(start), int(end)
            if start < 0 or end <= start:
                raise ValueError("Invalid BED interval")
            if intervals[chrom] and start < intervals[chrom][-1][1]:
                raise ValueError("BED intervals must be sorted and non-overlapping per contig")
            intervals[chrom].append((start, end))
    return dict(intervals)


def intersect_masks(left, right):
    shared = {}
    for chrom in sorted(left.keys() & right.keys()):
        a, b = left[chrom], right[chrom]
        i = j = 0
        hits = []
        while i < len(a) and j < len(b):
            start, end = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
            if start < end:
                hits.append((start, end))
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
        if hits:
            shared[chrom] = hits
    return shared


def read_snps(path, mask, excluded_contigs=()):
    result = {}
    starts = {chrom: [start for start, _ in intervals] for chrom, intervals in mask.items()}
    opener = gzip.open if str(path).endswith(".gz") else open
    total_alternate_snps = 0
    with opener(path, "rt") as stream:
        for line in stream:
            if line.startswith("#CHROM") and len(line.rstrip().split("\t")) != 10:
                raise ValueError("Expected a single-sample VCF")
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 10:
                raise ValueError("Expected a single-sample VCF")
            chrom, pos, _, ref, alt, _, filt, _, fmt, sample = fields
            if chrom in excluded_contigs:
                continue
            if filt not in {"PASS", "."} or len(ref) != 1 or len(alt) != 1:
                continue
            if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                continue
            values = dict(zip(fmt.split(":"), sample.split(":")))
            gt = values.get("GT", ".").replace("|", "/").split("/")
            if "." in gt or "1" not in gt:
                continue
            if any(allele not in {"0", "1"} for allele in gt):
                raise ValueError("Normalize and split multiallelic records before comparison")
            total_alternate_snps += 1
            position = int(pos)
            if chrom not in starts:
                continue
            index = bisect.bisect_right(starts[chrom], position - 1) - 1
            if index >= 0 and mask[chrom][index][0] <= position - 1 < mask[chrom][index][1]:
                key = (chrom, position, ref, alt)
                genotype = tuple(sorted(int(x) for x in gt))
                if key in result and result[key] != genotype:
                    raise ValueError("Conflicting duplicate SNP genotypes")
                result[key] = genotype
    return result, total_alternate_snps


def compare(illumina, ont, mask, excluded_contigs=()):
    mask = {chrom: intervals for chrom, intervals in mask.items() if chrom not in excluded_contigs}
    a, total_a = read_snps(illumina, mask, excluded_contigs)
    b, total_b = read_snps(ont, mask, excluded_contigs)
    shared = a.keys() & b.keys()
    union = a.keys() | b.keys()
    return {
        "excluded_contigs": sorted(excluded_contigs),
        "shared_depth_eligible_bases": sum(end - start for ranges in mask.values() for start, end in ranges),
        "illumina_alternate_snps_all_regions": total_a,
        "ont_alternate_snps_all_regions": total_b,
        "illumina_snps_in_shared_mask": len(a), "ont_snps_in_shared_mask": len(b),
        "shared_snps": len(shared), "illumina_only_snps": len(a.keys() - b.keys()),
        "ont_only_snps": len(b.keys() - a.keys()),
        "snp_jaccard": len(shared) / len(union) if union else None,
        "genotype_agreement_on_shared_snps": sum(a[k] == b[k] for k in shared) / len(shared) if shared else None,
        "interpretation": "Set overlap of non-reference SNPs, not precision/recall or genome-wide genotype concordance. A depth mask is not a validated callable-region mask."
    }, a, b


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("illumina", "ont", "illumina-bed", "ont-bed", "json", "tsv", "shared-bed"):
        p.add_argument("--" + name, required=True)
    p.add_argument("--exclude-contigs", nargs="*", default=[])
    args = p.parse_args()
    mask = intersect_masks(read_bed(args.illumina_bed), read_bed(args.ont_bed))
    mask = {chrom: intervals for chrom, intervals in mask.items() if chrom not in args.exclude_contigs}
    result, left, right = compare(args.illumina, args.ont, mask, args.exclude_contigs)
    with open(args.json, "w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    with open(args.shared_bed, "w") as stream:
        for chrom, ranges in mask.items():
            for start, end in ranges:
                stream.write(f"{chrom}\t{start}\t{end}\n")
    with open(args.tsv, "w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["chrom", "pos_1based", "ref", "alt", "category", "illumina_gt", "ont_gt"])
        for key in sorted(left.keys() | right.keys()):
            category = "shared" if key in left and key in right else "illumina_only" if key in left else "ont_only"
            writer.writerow([*key, category, "/".join(map(str, left[key])) if key in left else "",
                             "/".join(map(str, right[key])) if key in right else ""])
