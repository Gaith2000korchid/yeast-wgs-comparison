"""Controlled ONT base-quality sensitivity; reuse the completed primary BAMs.

Only the ONT base-quality cutoff changes. Illumina calls and its depth mask
remain exactly those from the primary analysis. No alignment is repeated.
"""
from pathlib import Path
import sys
configfile: "config/config.yaml"
sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from validate_inputs import validate_config
ROWS = validate_config(config)
SAMPLES = sorted({sample for sample, platform in ROWS})
OUT = config.get("output_dir", "results")
LOG = config.get("log_dir", "logs")
BQ = config.get("sensitivity_ont_base_quality", 7)
if type(BQ) is not int or BQ < 0 or BQ == config["min_base_quality"]:
    raise ValueError("sensitivity_ont_base_quality must be a distinct nonnegative integer")
DEST = f"{OUT}/sensitivity/ont_bq{BQ}"
REF = config["reference"]
SCRIPTS = Path(workflow.basedir) / "scripts"
wildcard_constraints:
    sample="[A-Za-z0-9][A-Za-z0-9_.-]*"

rule all:
    input:
        expand(f"{DEST}/{{sample}}.comparison.json", sample=SAMPLES),
        expand(f"{DEST}/{{sample}}.ont.pass.vcf.gz.csi", sample=SAMPLES)

rule ont_sensitivity_depth:
    input:
        bam=f"{OUT}/bam/{{sample}}.ont.bam",
        bai=f"{OUT}/bam/{{sample}}.ont.bam.bai",
        script=str(SCRIPTS / "coverage.py")
    output:
        bed=f"{DEST}/{{sample}}.ont.depth_mask.bed",
        json=f"{DEST}/{{sample}}.ont.coverage.json"
    params: mq=config["min_mapping_quality"], low=config["min_depth"], high=config["max_depth"]
    log: f"{LOG}/sensitivity/{{sample}}.ont_bq{BQ}.depth.log"
    conda: "../environment.yml"
    shell:
        "(samtools depth -aa -s -q {BQ} -Q {params.mq} -G 3844 {input.bam:q} | "
        "python {input.script:q} --min-depth {params.low} --max-depth {params.high} --bed {output.bed:q} --json {output.json:q}) 2> {log:q}"

rule ont_sensitivity_call:
    input:
        ref=REF, fai=f"{REF}.fai",
        bam=f"{OUT}/bam/{{sample}}.ont.bam", bai=f"{OUT}/bam/{{sample}}.ont.bam.bai"
    output:
        raw=f"{DEST}/{{sample}}.ont.raw.vcf.gz",
        vcf=f"{DEST}/{{sample}}.ont.pass.vcf.gz",
        index=f"{DEST}/{{sample}}.ont.pass.vcf.gz.csi",
        stats=f"{DEST}/{{sample}}.ont.bcftools.stats.txt"
    params:
        ploidy=lambda w: ROWS[w.sample, "ont"]["ploidy"],
        mq=config["min_mapping_quality"], cap=config["pileup_max_depth"],
        expr=f'QUAL="." || QUAL<{config["min_variant_quality"]} || FMT/DP="." || FMT/DP<{config["min_depth"]} || FMT/DP>{config["max_depth"]}'
    log: f"{LOG}/sensitivity/{{sample}}.ont_bq{BQ}.call.log"
    conda: "../environment.yml"
    shell:
        """
        (bcftools mpileup -B --skip-indels -q {params.mq} -Q {BQ} -d {params.cap} --ns 3844 \
          -a FORMAT/DP,FORMAT/AD -Ou -f {input.ref:q} {input.bam:q} | \
          bcftools call --ploidy {params.ploidy} -mv -Oz -o {output.raw:q}) 2> {log:q}
        (bcftools norm -f {input.ref:q} -c e -m -any -Ou {output.raw:q} | \
          bcftools view -e {params.expr:q} -Oz -o {output.vcf:q}) 2>> {log:q}
        bcftools index -f {output.vcf:q} 2>> {log:q}
        bcftools stats {output.vcf:q} > {output.stats:q} 2>> {log:q}
        """

rule compare_sensitivity:
    input:
        illumina=f"{OUT}/variants/{{sample}}.illumina.pass.vcf.gz",
        ont=f"{DEST}/{{sample}}.ont.pass.vcf.gz",
        illumina_bed=f"{OUT}/qc/{{sample}}/illumina.depth_mask.bed",
        ont_bed=f"{DEST}/{{sample}}.ont.depth_mask.bed",
        script=str(SCRIPTS / "compare_snps.py")
    output:
        json=f"{DEST}/{{sample}}.comparison.json",
        tsv=f"{DEST}/{{sample}}.snps.tsv",
        bed=f"{DEST}/{{sample}}.shared_depth_mask.bed"
    params: exclude=config.get("comparison_exclude_contigs", [])
    log: f"{LOG}/sensitivity/{{sample}}.ont_bq{BQ}.compare.log"
    conda: "../environment.yml"
    shell:
        "python {input.script:q} --illumina {input.illumina:q} --ont {input.ont:q} "
        "--illumina-bed {input.illumina_bed:q} --ont-bed {input.ont_bed:q} "
        "--json {output.json:q} --tsv {output.tsv:q} --shared-bed {output.bed:q} --exclude-contigs {params.exclude:q} 2> {log:q}"
