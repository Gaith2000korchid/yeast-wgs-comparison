"""Validate matched biological samples before constructing the workflow DAG."""
import csv
import re
from pathlib import Path


def load_samples(path):
    with open(path, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["sample", "platform", "read1", "read2", "ploidy"]:
            raise ValueError("Samples TSV must have: sample, platform, read1, read2, ploidy (in order)")
        rows = {}
        for row in reader:
            sample, platform = row["sample"], row["platform"]
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
                raise ValueError(f"Invalid sample identifier: {sample!r}")
            if platform not in {"illumina", "ont"} or row["ploidy"] not in {"1", "2"}:
                raise ValueError(f"Unsupported platform/ploidy for {sample}")
            if not row["read1"] or (platform == "illumina" and not row["read2"]):
                raise ValueError(f"Missing reads for {sample}/{platform}")
            if platform == "ont" and row["read2"]:
                raise ValueError("ONT reads must be single-end")
            if row["read1"] == row["read2"]:
                raise ValueError("Illumina R1 and R2 must be different files")
            key = (sample, platform)
            if key in rows:
                raise ValueError(f"Duplicate sample/platform row: {key}")
            rows[key] = row
        if not rows:
            raise ValueError("Empty samples TSV")
        for sample in {s for s, _ in rows}:
            if any((sample, p) not in rows for p in ("illumina", "ont")):
                raise ValueError(f"Both platforms are required for {sample}")
            if rows[sample, "illumina"]["ploidy"] != rows[sample, "ont"]["ploidy"]:
                raise ValueError(f"Ploidy mismatch for {sample}")
    return rows


def validate_config(config):
    required = {"reference", "samples", "dataset_label", "minimap2_preset",
                "min_mapping_quality", "min_base_quality", "min_depth", "max_depth",
                "min_variant_quality", "pileup_max_depth"}
    if required - config.keys():
        raise ValueError(f"Missing config keys: {sorted(required - config.keys())}")
    if config["minimap2_preset"] not in {"map-ont", "lr:hq"}:
        raise ValueError("Use map-ont for noisy ONT, lr:hq for high-quality ONT")
    for key in ("min_mapping_quality", "min_base_quality", "min_depth", "max_depth",
                "min_variant_quality", "pileup_max_depth"):
        if type(config[key]) is not int or config[key] < 0:
            raise ValueError(f"{key} must be a nonnegative integer")
    if not 0 < config["min_depth"] <= config["max_depth"] < config["pileup_max_depth"]:
        raise ValueError("Require 0 < min_depth <= max_depth < pileup_max_depth")
    reference = Path(config["reference"])
    if not reference.is_file():
        raise ValueError(f"Missing reference {reference}; generate demo or prepare real data first")
    rows = load_samples(config["samples"])
    for row in rows.values():
        for field in ("read1", "read2"):
            if row[field] and not Path(row[field]).is_file():
                raise ValueError(f"Missing input: {row[field]}")
    return rows
