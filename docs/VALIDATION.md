# Workflow validation

## Executed locally

- Fourteen Python unittest cases passed after adding real-preparation and contig-exclusion checks.
- Snakemake 9.13.4 constructed the complete 20-job demonstration DAG.
- The full synthetic analysis (both read branches, QC, BAMs, VCFs, comparison and reports) completed in the pinned environment.
- Integration assertion: all three planted SNPs recovered by both branches.
- Observed toy result: 3 shared SNPs, 0 Illumina-only, 0 ONT-only; shared depth-eligible region = 19,099 of 20,000 synthetic bases. These are synthetic software-test results only.

The local execution environment does not expose child process IDs to psutil in the usual way. Snakemake resource benchmarking failed with `psutil.NoSuchProcess`; local execution therefore used an otherwise identical temporary Snakefile with the benchmark directives omitted. Scientific commands and dependencies were unchanged. Benchmark collection remains enabled in the committed workflow and is covered by the GitHub-hosted synthetic integration job. Local validation alone does not establish that resource benchmarking works in this environment.

## Automated GitHub validation

The workflow runs unit tests and the committed synthetic workflow on an Ubuntu runner. The Actions status is authoritative for that run. A successful synthetic job proves execution and recovery of toy variants, not real-data biological correctness.

## Still required

- Real P11 execution and provenance verification completed; continue scientific validation beyond the baseline.
- Review coverage/filter choices and discordant loci in IGV.
- Add a suitable ONT-specific caller after reviewing chemistry and basecaller compatibility.
- Add a resolved environment lockfile and validated repeat/mappability exclusions if the project becomes a benchmark.

## Real-data preparation and QC regression

The preparation tests require synchronized Illumina mates, deterministic pair selection, valid FASTQ records and archive checksums. Historical ONT tests retain one consensus 2D record instead of correlated template/complement reads and reject duplicate consensus molecule identifiers. Comparison tests also cover excluded mitochondrial contigs.

The first real NanoPlot run exposed an upstream incompatibility between NanoPlot 1.43 and Plotly 6/7 when more than 10,000 reads are plotted ([upstream issue #400](https://github.com/wdecoster/NanoPlot/issues/400)). The environment explicitly pins Plotly 5.24.1 and pandas 2.2.3. CI includes a separate 12,001-read synthetic QC regression because the small workflow fixture could not expose that failure.


## Completed real execution

The primary P11 workflow completed successfully on matched public reads; its BAM integrity checks, filtered VCF indexes, QC and reports were produced. The separate four-job sensitivity workflow reused exactly the same ONT BAM and Illumina VCF/mask while changing only ONT base quality from 13 to 7. It also completed successfully. [Measured results and execution evidence](results/P11/README.md) retain both domains.

A clean pinned environment was used after an in-place local dependency update retained inconsistent package files. Execution resumed from completed scientific outputs after the NanoPlot failures. Successful completion records versions/configuration; failed QC attempts are acknowledged in the run record. No biological result from a failed or partial QC output is published.

The sensitivity integration test recovers identical toy comparison metrics at BQ13 and BQ7 on synthetic high-quality reads. Unit and integration checks validate computation, not clinical/biological accuracy. The real analysis has no truth set, validated repeat/mappability mask or manual IGV confirmation.
