# Yeast WGS comparison

Reproducible **Snakemake** workflow for matched **Illumina paired-end** and **Oxford Nanopore** reads from *Saccharomyces cerevisiae*: read QC, alignment, coverage, small-variant calling and alternate-SNP overlap.

**Status: initial implementation.** The included demonstration uses deterministic synthetic DNA, not a yeast genome. Biological conclusions require a verified matched public dataset. No real-data results are claimed here.

## What this project does

| Branch | Steps | Output |
|---|---|---|
| Illumina | FastQC → fastp → BWA-MEM → samtools fixmate/markdup → bcftools | QC, indexed BAM, normalized filtered SNP/indel VCF |
| Nanopore | NanoPlot → minimap2 → samtools → bcftools SNP baseline | QC, indexed BAM, normalized filtered SNP VCF |
| Comparison | Shared depth-eligible regions → exact alternate-SNP intersection | JSON metrics, SNP TSV, shared BED, HTML report |
| Reproducibility | Pinned tool versions, tool/version capture, reference SHA256, logs, benchmarks, CI | Traceable execution |

The ONT caller is an **exploratory bcftools baseline**, not a validated ONT calling solution. Indels and structural variants are not compared. A technology-appropriate caller can be added after chemistry/basecaller metadata are established.

## Run the synthetic demonstration

Use Linux, macOS, or Ubuntu under WSL2 on Windows. Install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) first.

```bash
git clone https://github.com/Gaith2000korchid/yeast-wgs-comparison.git
cd yeast-wgs-comparison
micromamba create -y -f environment.yml
micromamba activate yeast-wgs

python workflow/scripts/generate_demo.py
python -m unittest discover -s tests -v
snakemake -s workflow/Snakefile --cores 4 --dry-run
snakemake -s workflow/Snakefile --cores 4 --printshellcmds
python tests/check_demo.py
```

Open `results/report.html` and `results/multiqc/multiqc_report.html`. The demo plants three known SNPs and checks that both branches recover them. It tests workflow execution, not performance on real sequencing data. It has no realistic ONT indel or homopolymer error model.

The environment pins direct package versions; a resolved conda lockfile is not yet provided. CI also tests the same synthetic workflow from a clean checkout.

## Use public yeast reads

1. Follow [dataset selection](docs/DATASETS.md). Verify the **same isolate/sample** across technologies, organism, WGS strategy, library layout, chemistry/basecaller, ploidy and file checksums.
2. Copy `config/real.example.yaml` to `config/real.yaml`, replace the reference and sample sheet paths, and fill a copy of `config/samples.real.example.tsv`.
3. Keep the same `sample` identifier on both rows only when the biological match is verified. Supply an uncompressed FASTA reference and `.fastq` / `.fastq.gz` reads. For multiple lanes, concatenate lanes for each mate and retain the lane/run provenance separately.
4. Set ploidy to the documented value (1 or 2). The example's value of 1 is a placeholder, not an assertion about a strain. Aneuploid/polyploid genomes need a separate design.
5. Select `map-ont` for noisy ONT data or `lr:hq` for high-quality ONT data, following minimap2 documentation. Adjust depth thresholds to measured coverage.

```bash
snakemake -s workflow/Snakefile --configfile config/real.yaml --cores 4 --dry-run
snakemake -s workflow/Snakefile --configfile config/real.yaml --cores 4 --printshellcmds
```

## Read the results correctly

- **SNP Jaccard** = shared alternate SNPs / union of alternate SNPs, within the shared depth mask. Empty denominators produce `null`, not 100% agreement.
- **Genotype agreement** is computed only at shared alternate SNPs. It is not genome-wide genotype concordance. Missing and reference-only genotypes are excluded.
- A depth-eligible mask restricts the comparison to positions with sufficient, bounded high-quality coverage in both datasets. It does not prove that all those positions are callable: repeat/mappability and other quality exclusions are still needed.
- An Illumina-only or ONT-only SNP is **discordant**, not automatically true or false. No precision/recall claim is valid without a suitable truth set.
- Normalization/splitting precedes comparison; exact allele matching is used for SNPs. Indel representation equivalence is outside this implementation.

See [methodology and limitations](docs/METHODS.md) and the [French walkthrough](docs/DEMARRER_FR.md).

## Repository layout

| Path | Purpose |
|---|---|
| `workflow/Snakefile` | Both read branches and report dependencies |
| `workflow/scripts/` | Validation, demo generation, streaming coverage, SNP comparison, reporting |
| `config/` | Demonstration configuration and real-data templates |
| `tests/` | Coordinate, comparison, input-validation and integration checks |
| `docs/` | Methods, data provenance and learning sequence |
| `.github/workflows/ci.yml` | Unit tests and synthetic end-to-end execution |

Large reads, BAMs and generated outputs are ignored by Git. Publish small results with their provenance only after running and checking the real analysis.

## Learning sources and attribution

This workflow was written for this repository, informed by these resources:

- [Official Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) and [example data](https://github.com/snakemake/snakemake-tutorial-data): workflow dependencies and variant-calling example.
- [Chiara Batini, EMBL-EBI training materials (November 2025)](https://github.com/cbatini/training_materials/blob/main/EBI_NGS_Nov2025/days1_2_mapping_variant_calling_handbook_Nov2025.md): yeast Illumina QC and variant-calling practical.
- [BWA](https://github.com/lh3/bwa), [minimap2](https://github.com/lh3/minimap2), [SAMtools](https://www.htslib.org/doc/samtools.html), [BCFtools](https://www.htslib.org/doc/bcftools.html), [fastp](https://github.com/OpenGene/fastp), [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [NanoPlot](https://github.com/wdecoster/NanoPlot), [MultiQC](https://docs.seqera.io/multiqc): tool documentation.

MIT license applies to this repository's code. Third-party tools and datasets retain their own licenses. Implementation developed with AI assistance; results require execution and scientific review.
