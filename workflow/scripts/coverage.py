"""Stream samtools depth -aa to a depth mask and summary (0-based half-open BED)."""
import argparse
import json
import sys
from collections import defaultdict


def summarize(stream, bed, minimum, maximum):
    total = covered = eligible = depth_sum = 0
    by_chrom = defaultdict(lambda: {"bases": 0, "covered_bases": 0, "depth_eligible_bases": 0, "depth_sum": 0})
    start = end = chrom = None
    previous_chrom, previous_pos = None, 0
    completed = set()
    for line in stream:
        current, pos, depth = line.rstrip().split("\t")
        pos, depth = int(pos), int(depth)
        if pos < 1 or depth < 0:
            raise ValueError("Depth positions must be positive and depths nonnegative")
        if current != previous_chrom:
            if current in completed:
                raise ValueError("Depth contigs must not be interleaved")
            if previous_chrom is not None:
                completed.add(previous_chrom)
            previous_pos = 0
        if pos != previous_pos + 1:
            raise ValueError("Expected every position from samtools depth -aa, in order")
        previous_chrom, previous_pos = current, pos
        good = minimum <= depth <= maximum
        total += 1
        covered += depth > 0
        eligible += good
        depth_sum += depth
        item = by_chrom[current]
        item["bases"] += 1
        item["covered_bases"] += depth > 0
        item["depth_eligible_bases"] += good
        item["depth_sum"] += depth
        if good and chrom == current and end == pos - 1:
            end = pos
        else:
            if start is not None:
                bed.write(f"{chrom}\t{start}\t{end}\n")
            chrom, start, end = (current, pos - 1, pos) if good else (None, None, None)
    if start is not None:
        bed.write(f"{chrom}\t{start}\t{end}\n")
    if total == 0:
        raise ValueError("Empty depth stream")
    for item in by_chrom.values():
        item["mean_depth"] = item.pop("depth_sum") / item["bases"]
    return {"reference_bases": total, "covered_bases": covered,
            "depth_eligible_bases": eligible, "mean_depth": depth_sum / total,
            "covered_fraction": covered / total, "depth_eligible_fraction": eligible / total,
            "min_depth": minimum, "max_depth": maximum, "by_chromosome": dict(by_chrom)}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bed", required=True)
    p.add_argument("--json", required=True)
    p.add_argument("--min-depth", required=True, type=int)
    p.add_argument("--max-depth", required=True, type=int)
    a = p.parse_args()
    with open(a.bed, "w") as bed:
        result = summarize(sys.stdin, bed, a.min_depth, a.max_depth)
    with open(a.json, "w") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
