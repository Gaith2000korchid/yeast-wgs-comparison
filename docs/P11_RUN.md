# Reproduce the real P11 analysis

The biological identity is **CIC / Ponton11**, ENA sample **SAMEA3895683 / ERS1082817**, study **ERP016443**. [Recorded evidence](metadata/P11.provenance.json) links both technologies to this sample and the publication run table. The same named strain is documented as diploid, euploid and heterozygous in Peter et al. (2018), Supplementary Table S1. This supports a diploid nuclear calling model; it is not a new physical measurement of the archived preparation.

## Inputs and preparation

- Illumina HiSeq 2000: ERR1527933 paired FASTQs, approximately 2.90 GB compressed.
- Historical MinION: ERR1539063–ERR1539067, approximately 0.54 GB compressed. ENA labels the experiment R7; the study mostly used R7.3. Per-run exact chemistry/basecaller versions were not established.
- Reference: S288C R64-1-1, Ensembl release 113, 16 nuclear chromosomes plus Mito. Its pinned URL and both compressed/uncompressed SHA256 digests are in `config/P11.reference.json`.
- The reference digests were computed after the initial HTTPS download and gzip CRC check; they are reproducibility fingerprints, not independent publisher-authenticated SHA256 values.

Download requires about 3.44 GB, plus reference. Reserve at least 15 GB free for preparation, reads, alignments, temporary sorting and reports. Downloads are excluded from Git.

From the repository root, in the `yeast-wgs` environment:

```bash
python workflow/scripts/download_reference.py
python workflow/scripts/download_reads.py --manifest config/P11.downloads.json --workers 2
python workflow/scripts/prepare_real_reads.py --fraction 0.25 --seed 20261001
snakemake -s workflow/Snakefile --configfile config/P11.yaml --cores 4 --dry-run
snakemake -s workflow/Snakefile --configfile config/P11.yaml --cores 4 --printshellcmds
# Controlled sensitivity: only ONT BQ changes, from 13 to 7.
snakemake -s workflow/quality_sensitivity.smk --configfile config/P11.yaml --cores 2 --printshellcmds
```

The downloader validates archive byte counts and MD5 before atomically completing each file. Reruns verify completed files and skip them. Failed partial downloads restart; there is no range-resume implementation.

The preparer revalidates all archived files, parses complete FASTQ records and requires matching Illumina read identifiers and equal mate counts. It selects each pair if the first 64 bits of `SHA256("20261001:" + read_name)` fall below 25% of the 64-bit range. Both mates share the decision. This is a random-like full-genome subset, not the first quarter of the file. Gzip output uses fixed timestamps. Exact selected read/base counts and output SHA256 digests are recorded in `data/real/prepared/preparation.json`.

## Historical ONT read selection

These archive FASTQs contain template, complement and consensus (`twodirections`) records from the same sequencing molecule. Combining every record would inflate coverage with correlated observations. The preparation therefore selects **only consensus 2D records, both historical pass and fail**. It records all six class counts per run and rejects repeated consensus molecule identifiers across runs. It does not select only the smallest run or silently replace missing consensus reads with template reads.

This yields about 168 Mb of 2D sequence, not the approximately 598 Mb in all archived records. Historical pass/fail is a basecaller read classification, not a SAM QC-failure flag. Both classes remain for this documented baseline; mapping/base quality filters act downstream. Modern ONT-trained callers must not be assumed compatible with these old reads.

## Outputs and interpretation

The dedicated configuration writes `results/P11/`, `logs/P11/` and `benchmarks/P11/`, keeping demo outputs separate. Open `results/P11/report.html` and its MultiQC link. Main SNP comparison excludes `Mito`: nuclear diploidy is not an appropriate model for mitochondrial inheritance/copy number. BAMs and coverage QC retain all reference contigs; VCFs also retain mitochondrial calls, which are not biologically interpreted here.

Default thresholds are MAPQ ≥20, base quality ≥13, depth 8–200 and variant QUAL ≥30. The comparison domain is the intersection of depth-eligible nuclear bases, with each VCF additionally filtered by its own FORMAT/DP. These are exploratory starting thresholds, not validated callable-region criteria.

The two inputs have unequal coverage and different error profiles. SNP overlap is conditional on this mask and this bcftools baseline; it is not accuracy, sensitivity or a platform ranking. Repeat/mappability masks, an independent truth set and manual discordance review would be necessary for a defensible accuracy benchmark.

## Local execution limitation

The managed execution environment used for this project cannot expose child PIDs to psutil, causing Snakemake resource benchmarking to fail. Real and synthetic local runs use a temporary otherwise identical Snakefile with `benchmark:` lines omitted; no scientific command or dependency is changed. The committed workflow retains benchmarks, and GitHub's synthetic CI executes them on an ordinary Ubuntu runner. No local CPU/memory benchmark values are claimed.

## Export compact evidence

After completion, `python workflow/scripts/summarize_real_run.py` exports measured QC, coverage, SNP overlap and review loci to `docs/results/P11/`. This optional export step requires `matplotlib==3.10.8`; install it separately if absent. Raw reads, BAMs, VCFs and the full MultiQC HTML remain local and can be regenerated with the commands above. The GitHub archive contains compact metrics and provenance.

[Completed measured results](results/P11/README.md) document that the strict Q13 domain is very small and show the broader Q7 experiment. Their domains differ; their Jaccard values cannot be used as an accuracy ranking.
