# Real yeast WGS results: CIC / Ponton11 (P11)

Executed on **2026-10-01**, using matched ENA sample **SAMEA3895683**, reference **S288C R64-1-1 / Ensembl 113** and a documented **diploid** nuclear model. [Input identity and ploidy evidence](../../metadata/P11.provenance.json) · [Reproduce the run](../../P11_RUN.md) · [Machine-readable metrics](metrics.json).

The primary workflow and a controlled ONT base-quality sensitivity workflow completed successfully. These are measured results from real sequencing reads, separate from the synthetic CI demonstration.

## Read and alignment QC

| Measurement | Illumina | Historical ONT |
|---|---:|---:|
| Prepared input | 2,105,894 paired fragments | 29,914 consensus 2D reads |
| Input bases | 959,448,222 | 167,666,033 |
| Primary reads aligned | 97.40% | 74.64% |
| Primary reads marked duplicate | 0.25% | Not marked for ONT |
| Mean nuclear depth, MAPQ≥20 / BQ≥13 | 45.93× | 1.72× |
| Mean nuclear depth, ONT BQ≥7 sensitivity | Unchanged | 8.43× |

Illumina was selected using a deterministic 25% **paired** subset, seed 20261001. fastp retained all selected reads and trimmed approximately 0.044% of bases; post-processing Q30 base fraction was 90.94%. FastQC reported a **per-sequence GC-content failure** for both mates before and after fastp, and a length-distribution warning. Per-base quality and adapter content passed. These flags are retained in [module outcomes](fastqc_checks.tsv); they are not proof of contamination or a reason to discard the sample without further investigation. Taxonomic contamination analysis was not performed.

Historical ONT mean read quality was **Q8.0**, median Q8.8, read N50 **7,754 bp** ([NanoStats](NanoStats.txt)). Original archive reads mixed template/complement/consensus records from the same molecule; only consensus 2D pass **and** fail records were retained. This reduced approximately 598 Mb of correlated archive observations to 167.7 Mb of consensus sequence. Original basecaller pass/fail is distinct from SAM QC-failure flags. These reads are not modern high-accuracy ONT data.

![Nuclear filtered coverage and depth eligibility](coverage.svg)

The figure retains zero-coverage bases in chromosome means. Mito is omitted from nuclear metrics; its QC remains in [per-contig coverage](coverage_by_chromosome.tsv). Its copy number and inheritance require a different calling model. Coverage alone does not establish aneuploidy.

## SNP overlap and the quality-threshold experiment

Illumina calls and depth mask are **identical** in the two comparisons. Only the ONT minimum base-quality threshold changes from 13 to 7. MAPQ≥20, depth 8–200, variant QUAL≥30, reference, BAM, diploid model and all other parameters remain fixed. The experiment measures sensitivity to a filter choice; Q7 is not asserted to be an optimal or validated threshold.

| Metric | Primary: ONT BQ13 | Sensitivity: ONT BQ7 |
|---|---:|---:|
| Shared depth-eligible nuclear bases | 96,235 | 6,435,340 |
| Fraction of nuclear reference | 0.80% | 53.31% |
| Illumina alternate SNPs in shared mask | 296 | 38,412 |
| ONT alternate SNPs in shared mask | 160 | 5,471 |
| Shared alternate SNPs | 156 | 5,444 |
| Illumina-only alternate SNPs | 140 | 32,968 |
| ONT-only alternate SNPs | 4 | 27 |
| Alternate-SNP Jaccard | 52.00% | 14.16% |
| Genotype agreement at shared alternate SNPs only | 98.08% | 98.68% |

**Neither genotype-agreement percentage is genome-wide concordance or accuracy.** At BQ7 the caller produces only **8 heterozygous nuclear SNP records** out of 5,545 alternate SNP records across all retained nuclear regions; Illumina produces 35,684 heterozygous records out of 77,276. Shared ONT calls therefore strongly select a subset of largely homozygous alternate sites. The many Illumina-only sites remain relevant even though shared-site genotype agreement is high.

The two masks cover different genome domains, so comparing their Jaccard values as a platform ranking or as improvement/worsening in accuracy is invalid. The strict BQ13 result covers only 0.80% of nuclear bases. The broader BQ7 experiment still has substantial coverage and alternate-call imbalance. Old ONT errors, unequal effective coverage, mapping, filters and the exploratory bcftools caller can all contribute; this analysis does not separate their causal contributions or label discordant variants as false positives/negatives.

The same strain's diploidy is documented in an independent collection table. Illumina heterozygous allele balance concentrates near 0.5 ([histogram counts](metrics.json)), which is compatible with that model, but is not an independent ploidy measurement or a truth set. Published variant counts from other preparations/callers are not used as benchmark truth.

## Evidence, validation and next scientific checks

The archive contains exact preparation counts/checksums, [tool versions and resolved configuration](tool_versions.txt), [successful workflow record](run_record.json), [sensitivity record](sensitivity_run_record.json), [resolved Linux package URLs](environment-linux-64.explicit.txt), compact SAMtools/bcftools summaries and result fingerprints in [SHA256SUMS](SHA256SUMS). The Linux URL list can be used with `micromamba create -n yeast-wgs-replay -f environment-linux-64.explicit.txt`; it is platform-specific and has no independent package archive checksum manifest.

Local execution resumed after two NanoPlot QC failures and used a clean environment with pinned Plotly 5.24.1 / pandas 2.2.3. The committed environment and a >10,000-read regression are tested in CI. Local Snakemake resource benchmarks were omitted because managed child PIDs are not visible to psutil; no CPU/memory benchmark is claimed. Scientific commands/dependencies were unchanged by that omission.

Raw FASTQs, BAMs, full VCFs and full HTML QC outputs are excluded from Git and regenerated by the documented commands. The compact archive retains measured metrics, not an entire variant truth set. [Primary example loci](loci_for_review.tsv) and [BQ7 example loci](bq7.loci_for_review.tsv) are deterministic starting points for IGV review, **not manually validated variants**.

The [selected-locus follow-up](review/README.md) measures original-base support at 15 archived BQ7 loci and prepares a portable IGV session and review worksheet. A [targeted IGV Web review of six loci](review/IGV_WEB_REVIEW_FR.md) records reviewer-supplied coverage counts and two screenshot observations; nine selected loci remain unreviewed. No independent variant validation is claimed.

A defensible accuracy comparison needs repeat/mappability exclusions, suitable truth data, discordance review and a caller validated for the read chemistry. The present result demonstrates reproducible processing and exposes an important filter/caller limitation on legacy ONT data.
