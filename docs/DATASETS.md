# Public yeast dataset selection

## Selected real sample: CIC / Ponton11 (P11)

Istace et al. (2017), *De novo assembly and population genomic survey of natural yeast isolates with the Oxford Nanopore MinION sequencer*. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5466710/), [ENA study ERP016443](https://www.ebi.ac.uk/ena/browser/view/ERP016443).

The [saved ENA response](metadata/ERP016443.ena.tsv), retrieved 2026-10-01, contains paired Illumina and single-end ONT WGS FASTQ links for shared sample accessions. One relatively small candidate is **SAMEA3895683**, alias `sample_BCM_4932_YJS5845,P11`:

| Technology | Runs | Approximate compressed download |
|---|---|---|
| Illumina | ERR1527933 (R1/R2) | 2.90 GB |
| ONT | ERR1539063, ERR1539064, ERR1539065, ERR1539066, ERR1539067 | 0.54 GB |

Total: approximately **3.44 GB compressed**, with larger working storage required. Links, sizes and MD5 checksums are in the metadata snapshot. The shared sample accession and both experiment XML records verify the technology match. Istace et al. supplementary Table S5 identifies CIC as Ponton11, and Supplementary File 3 assigns all five ONT runs to CIC. Peter et al. (2018), Supplementary Table S1, identifies CIC / Ponton 11 as diploid, euploid and heterozygous. [Evidence](metadata/P11.provenance.json) and [run instructions](P11_RUN.md) record this chain. All seven files were downloaded and verified against ENA MD5/byte counts on 2026-10-01. These older ONT data are useful for learning but must not be treated as modern high-accuracy ONT reads.

Inspection of original FASTQ headers reveals mixed template/complement/consensus observations. The P11 preparer retains consensus 2D (`twodirections`) records only, both historical pass and fail, avoiding repeated observations from a molecule. Run and read-class counts are retained in preparation provenance. The smallest ONT run alone is not a substitute for sufficient coverage.

## Candidate with both read technologies

Giordano et al. (2017), *De novo yeast genome assemblies from MinION, PacBio and MiSeq platforms*, Scientific Reports 7, 3935. [DOI](https://doi.org/10.1038/s41598-017-03996-z), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5479803/).

The study reports four strains: S. cerevisiae S288C/SK1 and S. paradoxus N44/CBS432. The accession section identifies:

| Study | Data |
|---|---|
| [PRJEB19900](https://www.ebi.ac.uk/ena/browser/view/PRJEB19900) | ONT and MiSeq |
| [PRJEB7245](https://www.ebi.ac.uk/ena/browser/view/PRJEB7245) | PacBio |
| [S-BSST17](https://www.ebi.ac.uk/biostudies/studies/S-BSST17) | Generated assemblies |

This is a **candidate**, not a dataset already analyzed by this repository. These are older ONT reads and PacBio RSII data, not modern Q20 ONT or PacBio HiFi. Error profiles and caller/model compatibility must be considered. Select one S. cerevisiae sample with verified Illumina/ONT identity rather than comparing arbitrary runs labeled with the same species.

## Archive availability check (2026-10-01)

The [saved ENA metadata response](metadata/PRJEB19900.ena.tsv) was retrieved using the ENA `filereport` API for this study. It identifies Illumina run **ERR1938683** and several ONT runs for the same S288C sample accession **SAMEA4461732**. However, the ONT rows have **empty `fastq_ftp` fields**, while the Illumina row supplies two FASTQ links and MD5 checksums. The SK1 ONT rows also have empty FASTQ link fields.

Consequently, these S. cerevisiae data are not yet a ready-to-download paired FASTQ dataset for this workflow. Submitted-file availability, alternative mirrors or a newer matched study need checking. Do not invent FASTQ links by substituting an accession into another run's URL pattern. The S. paradoxus rows with FASTQ links do not meet this project's S. cerevisiae target.

Metadata query: `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB19900&result=read_run&fields=run_accession,sample_accession,sample_alias,scientific_name,instrument_platform,instrument_model,library_layout,library_strategy,fastq_ftp,fastq_md5,fastq_bytes&format=tsv`.

## Before a real run

Record the study/run/sample accessions, strain, organism, layout, WGS strategy, exact FASTQ links and archive MD5 checksums; confirm read file availability. Inspect the publication and supplementary material for sample identity, ploidy, ONT chemistry and basecaller history. Record reference accession/release and checksum, download date, any concatenation or subsampling commands and their seeds.

The example sample sheet intentionally contains placeholders. Do not infer ploidy from a familiar strain name alone. Do not substitute a different isolate to fill a missing technology.

For a small first analysis, use one verified isolate and the full reference. A later coverage-sensitivity experiment can use explicit, reproducible read subsampling. A one-chromosome read subset needs alignment-based read selection; truncating FASTQ records does not isolate a chromosome.
