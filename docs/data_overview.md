# Data overview
 
## Purpose
 
This document summarizes the datasets underlying the study *"Kin discrimination in Bacillus subtilis is predicted by average nucleotide identity but shaped by mobile genetic elements"* (Stare et al.) and clarifies which components are included in this repository and which are not.
 
---
 
## Study scope
 
- Bacillus subtilis genomes analyzed: 313 (40 closed PS-collection isolates, 6 additional newly sequenced wild-type strains, and 267 publicly available genomes)
- Pairwise swarming interactions scored: 850, across 78 strains
- Kin discrimination (KD) groups defined: 26
---
 
## Data included in this repository
 
This repository contains analysis code and documentation supporting the manuscript. It does not redistribute primary sequence data.
 
- Analysis and figure-generation scripts (S1–S18)
- Python and R package lists (`environment/`)
- Software, database, and environment versions (`docs/SOFTWARE_VERSIONS.md`)
The final curated supplementary tables are published with the article and are not duplicated here.
 
---
 
## Data not included in this repository
 
The following are not redistributed here, owing to their size and the computational infrastructure required to generate and process them:
 
- Raw sequencing reads
- Genome assemblies
- Raw outputs of bioinformatics tools
- Intermediate files and outputs
- Reference databases
- Manuscript figures, supplementary figures, and supplementary tables (available with the published article)
---
 
## Raw data provenance
 
Primary sequence data are available from public archives rather than this repository:
 
- Newly sequenced B. subtilis genomes — NCBI, under the accessions listed in the manuscript (see Methods and Supplementary Tables)
- Publicly available B. subtilis genomes (267) — NCBI, accessions listed in the Supplementary Information
- Prophage vOTU reference sequences — Štefanič et al. (2025), Zenodo 13834859
Software and database versions used to process these data are recorded in `docs/SOFTWARE_VERSIONS.md`.
 
---
 
## Repository data philosophy
 
This repository is structured as a documentation and transparency companion rather than a complete, re-executable pipeline. The goal is to provide:
 
- transparency of the analytical logic,
- documentation of software versions and the computational environment,
- and a record of manual curation decisions, rather than a full end-to-end re-executable workflow. Inputs are described in the manuscript Methods; but note that the scripts are not packaged for one-click execution.
 
---
 
## Notes on manual curation
 
Some inputs and tables were manually harmonized prior to analysis or inclusion. These operations included:
 
- reclassifying PHASTER prophage labels against vOTU references
- manually classifying BF (BLAST family) families from TAfinder output
- standardizing categorical labels and identifiers
- reformatting and harmonizing matrix structures for consistency
- renaming column headers
These curation steps are documented in the manuscript and its Supplementary Information.
