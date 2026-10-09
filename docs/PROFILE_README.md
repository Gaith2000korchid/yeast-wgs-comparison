# Gaith Korchid

Bioinformatics engineer · Biomedical data science · Scientific Python & R

I am a graduate of École d’Ingénieur Denis Diderot, Université Paris Cité, specialising in bioinformatics. I work at the intersection of biological data, quantitative methods and scientific software.

My research experience includes proteomic and clinical biomarker modelling at Brain&Mind, quantitative muscle MRI at Institut de Myologie, and genomics at Paris Brain Institute. The public projects below are separate portfolio studies: they document methods, contributions, results and reproducibility.

## Selected projects

### [Yeast WGS comparison](https://github.com/Gaith2000korchid/yeast-wgs-comparison)
A Snakemake workflow comparing matched Illumina and Nanopore yeast reads, from QC and alignment to filtered SNP comparison and targeted discordance review.

- **Evidence:** recorded real-data results, source fingerprints, unit tests and synthetic integration CI.
- **Finding:** comparison eligibility changes markedly with ONT base-quality thresholds; filtering explains selected discordances.
- **Tools:** Python, Snakemake, BWA, minimap2, samtools, bcftools.
- **Scope:** exploratory ONT calling without an independent truth set.

### [Reproducible cervical biopsy modelling](https://github.com/Gaith2000korchid/reproducible-cervical-biopsy-modeling)
A public-data modelling study with training-only preprocessing, nested group-aware validation, a dummy baseline and four sensitivity analyses.

- **Evidence:** out-of-fold predictions, tests, reproducible figures and a documented complete refit.
- **Finding:** modest primary discrimination and substantial inflation when current diagnostic variables are included. Numerical refit differences are documented.
- **Tools:** Python, scikit-learn, XGBoost, pandas, uv.
- **Scope:** internal methodological evaluation; no external clinical validation.

### [Biopython M-pro Atlas](https://github.com/Gaith2000korchid/biopython-mpro-atlas)
A CLI joining PROSITE, UniProtKB, PDB and SIFTS to inspect domain coverage and structural annotations.

- **Evidence:** dated source fingerprints, offline replay and tests under Python 3.10 and 3.12.
- **Finding:** an empty primary PROSITE PDB field does not mean structures are absent; identifier overlap needs residue-range checks.
- **Tools:** Python, Biopython, biological APIs.
- **Scope:** annotation and coordinate analysis, without docking or affinity prediction.

## Additional projects

| Project | What it demonstrates |
|---|---|
| [Viral genome comparator R](https://github.com/Gaith2000korchid/viral-genome-comparator-r) | Tested R package; a frozen 29-genome panel; alignment, identity/coverage, conservation and exploratory taxonomy analysis. |
| [De Bruijn genome assembler](https://github.com/Gaith2000korchid/debruijn-genome-assembler) | Attributed educational implementation; Python generators, weighted graphs, regression tests, Docker and FASTQ-to-FASTA CI. |

For genomics and R roles, start with Yeast WGS and Viral Genome Comparator. For Python development, the assembler provides complementary algorithmic evidence.

## Skills and interests

**Demonstrated in these repositories:** scientific Python and R, biological data processing, graph algorithms, predictive-model evaluation, workflow automation, testing, GitHub Actions, dependency management and execution provenance.

**Broader training and experience:** SQL, biomedical statistics, quantitative imaging and data platforms. I am also developing a data-engineering direction; the repositories above primarily demonstrate scientific computing and biological analysis.

I am interested in bioinformatics, biomedical data science, scientific software, research engineering and relevant R&D opportunities.

## Contact

[GitHub](https://github.com/Gaith2000korchid) · [Email](mailto:korchid2000gaith@gmail.com)
