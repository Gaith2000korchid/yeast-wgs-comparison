# Initial implementation validation

## Executed locally

- Nine Python unittest cases passed.
- Snakemake 9.13.4 constructed the complete 20-job demonstration DAG.
- The full synthetic analysis (both read branches, QC, BAMs, VCFs, comparison and reports) completed in the pinned environment.
- Integration assertion: all three planted SNPs recovered by both branches.
- Observed toy result: 3 shared SNPs, 0 Illumina-only, 0 ONT-only; shared depth-eligible region = 19,099 of 20,000 synthetic bases. These are synthetic software-test results only.

The local execution environment does not expose child process IDs to psutil in the usual way. Snakemake resource benchmarking failed with `psutil.NoSuchProcess`; local execution therefore used an otherwise identical temporary Snakefile with the benchmark directives omitted. Scientific commands and dependencies were unchanged. Benchmark collection remains enabled in the committed workflow and is covered by the GitHub-hosted synthetic integration job. Local validation alone does not establish that resource benchmarking works in this environment.

## Automated GitHub validation

The workflow runs unit tests and the committed synthetic workflow on an Ubuntu runner. The Actions status is authoritative for that run. A successful synthetic job proves execution and recovery of toy variants, not real-data biological correctness.

## Still required

- Verify and run a public matched yeast sample, including ploidy and raw-read provenance.
- Review coverage/filter choices and discordant loci in IGV.
- Add a suitable ONT-specific caller after reviewing chemistry and basecaller compatibility.
- Add a resolved environment lockfile and validated repeat/mappability exclusions if the project becomes a benchmark.
