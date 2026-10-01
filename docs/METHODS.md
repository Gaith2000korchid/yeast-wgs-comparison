# Methods and limits

## Input contract

A sample sheet contains exactly one Illumina paired-end row and one ONT row per biological sample, with matching ploidy (1 or 2). Identical identifiers are a declaration by the user, not automated evidence of sample identity. FASTQ basecalls are the starting point; POD5/FAST5 processing is outside scope.

The workflow rejects missing files, duplicated rows, incomplete platform pairs, malformed IDs and incompatible ploidy. Reference FASTA must be uncompressed and writable for its indexes. Rows may point to external local read files; paths are interpreted from the repository working directory.

## Read processing and mapping

Illumina QC is performed before/after fastp, with paired-end adapter detection. BWA-MEM maps the cleaned reads. Name sorting, fixmate, coordinate sorting and markdup mark duplicate fragments; duplicates are retained in the BAM but excluded from depth and calling. ONT basecalls are QC'd with NanoPlot and mapped using the configured minimap2 preset. ONT reads are not passed through Illumina-specific trimming or paired-end duplicate marking.

Read groups preserve the biological sample identifier and platform. SAMtools quickcheck checks basic BAM integrity; flagstat/stats and MultiQC expose mapping outcomes. A future iteration should include an explicit failure threshold for poor mapping, contamination/sample identity analysis and real-data review in IGV.

## Coverage mask

SAMtools depth `-aa` includes zero-coverage positions and unused reference contigs. Base and mapping quality thresholds are explicit. Flags 3844 (0xF04) exclude unmapped, secondary, QC-failed, duplicate and supplementary alignments. `-s` avoids counting overlapping mates twice. The streaming Python step summarizes mean depth (including zero bases), coverage fraction and depth-eligible fraction, with per-contig summaries.

Positions between the configured minimum and maximum depth are written to 0-based half-open BED intervals. The mask intersection is the comparison domain. The mask is a **depth mask**, not a benchmark-grade callable mask. Depth from SAMtools and FORMAT/DP from bcftools can differ because overlap/base handling differs; both restrictions are applied. Depth summaries do not establish copy-number changes or aneuploidy.

## Calling and filtering

BCFtools mpileup/call runs independently for each technology with the configured fixed ploidy. ONT uses `-B --skip-indels`, an exploratory SNP-only baseline that disables short-read-oriented BAQ correction. No ONT indel accuracy claim is made. A caller validated for the selected ONT data should be benchmarked separately before interpreting biological discordances.

Calls are normalized against the exact reference, multiallelic records split, and records with missing/low QUAL or missing/out-of-range FORMAT/DP excluded. Filter values are starting values, not universal validated cutoffs. SNP sets include only passing, biallelic A/C/G/T alternate alleles with complete non-reference genotypes. Genotypes are compared without phasing.

Illumina indels are retained in the filtered VCF and counted by bcftools stats, but are not cross-platform benchmarked. No structural-variant caller, annotation, phasing, joint cohort calling or ploidy inference is implemented.

## Metrics

Let A/B be the two alternate-SNP sets inside the shared depth mask. Jaccard = |A ∩ B| / |A ∪ B|. Shared-SNP genotype agreement is the fraction of identical unphased GTs among shared alleles only. Undefined fractions are null. The TSV retains shared and platform-only SNPs for later IGV review.

Neither metric is precision/recall. Even identical calls can share an error, and two different results can arise from coverage, mapping, repeats, genotype likelihoods or sequencing errors. Shared high-depth coverage is not a truth set.

## Reproducibility and validation

Direct package versions are pinned in environment.yml. The run captures resolved configuration, sample sheet, reference/sample-sheet SHA256 and tool versions. Raw-read checksums and ENA metadata must be kept in a real-data manifest; the first implementation does not hash large FASTQ files automatically. Logs/benchmarks are per sample/rule. The environment has no transitive dependency lockfile yet.

Tests cover BED boundaries, depth thresholds/zero bases, mask intersection, reference and missing genotypes, duplicate BED intervals, empty denominators, paired platform/ploidy validation and deterministic fixture generation. The integration check requires all three planted synthetic SNPs to be recovered by both branches; this is an execution check, not a realistic sequencing benchmark.
