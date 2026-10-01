"""Deterministic toy DNA/read fixtures. Not a yeast genome or an ONT simulator."""
import argparse
import gzip
import io
import json
import random
from pathlib import Path


def write_fastq(path, records):
    # Deterministic gzip bytes (no filename/timestamp in the header).
    with open(path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as zipped:
            with io.TextIOWrapper(zipped) as out:
                for name, seq, quality in records:
                    out.write(f"@{name}\n{seq}\n+\n{quality * len(seq)}\n")


def generate(directory):
    rng = random.Random(20261001)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    reference = "".join(rng.choice("ACGT") for _ in range(20000))
    mutations = {3001: "", 8001: "", 15001: ""}
    sample = list(reference)
    for pos in mutations:
        mutations[pos] = rng.choice([b for b in "ACGT" if b != reference[pos - 1]])
        sample[pos - 1] = mutations[pos]
    sample = "".join(sample)
    (directory / "reference.fa").write_text(
        ">synthetic_contig\n" + "\n".join(reference[i:i + 80] for i in range(0, len(reference), 80)) + "\n"
    )
    r1, r2 = [], []
    complement = str.maketrans("ACGT", "TGCA")
    for i in range(3000):
        start = rng.randrange(len(sample) - 450)
        r1.append((f"pair{i}/1", sample[start:start + 150], "I"))
        r2.append((f"pair{i}/2", sample[start + 300:start + 450].translate(complement)[::-1], "I"))
    long_reads = []
    for i in range(400):
        start = rng.randrange(len(sample) - 2000)
        seq = list(sample[start:start + 2000])
        for j, base in enumerate(seq):
            if rng.random() < 0.005:
                seq[j] = rng.choice([b for b in "ACGT" if b != base])
        long_reads.append((f"long{i}", "".join(seq), "5"))
    write_fastq(directory / "illumina_R1.fastq.gz", r1)
    write_fastq(directory / "illumina_R2.fastq.gz", r2)
    write_fastq(directory / "ont.fastq.gz", long_reads)
    (directory / "truth.json").write_text(json.dumps({
        "description": "Synthetic SNP fixture; no biological inference", "seed": 20261001,
        "snps": [{"chrom": "synthetic_contig", "pos": p, "ref": reference[p - 1], "alt": a}
                 for p, a in mutations.items()]
    }, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="data/demo")
    generate(parser.parse_args().out)
