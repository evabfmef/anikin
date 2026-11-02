# Data overview
 
## Purpose
 
This document summarizes the datasets underlying the study *"Genomic boundaries mark social boundaries in Bacillus subtilis"* (Stare et al.) and clarifies which components are included in this repository and which are not.
 
---
 
## Study scope
 
- Bacillus subtilis genomes analyzed: 313 (40 closed PS-collection isolates, 6 additional newly sequenced wild-type strains, and 267 publicly available genomes)
- Genomovar groups (≥99.5% ANI): 94
- Pairwise swarming interactions scored: 850, across 78 strains
- Kin discrimination (KD) groups: 25 candidate KD groups, refined to 26 final KD groups (67 strains)
- KD groups of the 39-strain PS panel: 14
---
 
## Data included in this repository
 
This repository contains analysis code and documentation supporting the manuscript. It does not redistribute primary sequence data.
 
| Folder | Contents |
| --- | --- |
| `scripts/` | Analysis and figure-generation scripts S1–S20 and the shared configuration loader `kd_config.py` (see `scripts/README.md`) |
| `config/` | Analysis settings: `config_S15.yml` (S15A, S15B) and `config_S16_S19.yml` (S17A–S17D, S19) |
| `metadata/` | Small strain-level tables read by the scripts (see below) |
| `environment/` | Software environments: `environment.yml` / `requirements.txt` (all Python scripts except S11B and S11C), `requirements_GEM.txt` (S11B, S11C) and `R_packages.txt` |
| `docs/` | This file, `software_versions.md` and `script_mapping.md` (inputs, outputs and dependencies of each script) |

### Metadata tables
 
| File | Contents | Used by |
| --- | --- | --- |
| `accessions.tsv` | All 313 genomes: strain name, NCBI assembly and sequence accession, genome file name, public/in-house flag | S17B, S17C (strain name ↔ genome file) |
| `kd_groups_39.tsv`, `kd_groups_67.tsv` | KD group membership of the 39- and 67-strain panels | S17A, S17D |
| `kd_groups_representatives_39.tsv`, `kd_groups_representatives_67.tsv` | Up to three representative strains per KD group used for variant validation | S17B |
 
KD group labels follow the manuscript (KD1–KD14 for the 39-strain panel; KD1–KD97 for the 67-strain panel). Labels are panel-specific: the same label can denote different groups in the two panels.
 
The final curated supplementary tables are published with the article and are not duplicated here. Scripts S6 and S14 read directly from the published Supplementary Tables workbook (S6: sheet 8; S14: sheets 5, 9 and 10).

---
 
## Data not included in this repository
 
The following are not redistributed here, owing to their size and the computational infrastructure required to generate and process them:
 
- Raw sequencing reads
- Genome assemblies and annotations
- Raw outputs of bioinformatics tools
- Intermediate files and outputs (feature matrices, interaction tables)
- Reference databases
- Manuscript figures, supplementary figures, and supplementary tables (available with the published article)

Input paths in the scripts and configuration files (e.g. `data/strain_interactions.xlsx`, `data/kd/inputs/`, `data/raw/gbk/`) refer to the local working layout used for the analysis. These files are not part of the repository; their content and origin are described in the manuscript Methods and in `docs/script_mapping.md`.

---
 
## Raw data provenance
 
Primary sequence data are available from public archives rather than this repository:
 
- Newly sequenced *B. subtilis* genomes: NCBI, under the accessions listed in the manuscript (see Methods and Supplementary Tables)
- Publicly available *B. subtilis* genomes (267): NCBI, accessions listed in `metadata/accessions.tsv` and in the Supplementary Information
- Prophage vOTU reference sequences: Štefanič et al. (2025), Zenodo 13834859
Software and database versions used to process these data are recorded in `docs/software_versions.md`.
 
---
 
## Repository data philosophy
 
This repository is structured as a documentation and transparency companion rather than a complete, re-executable pipeline. It provides:
 
- transparency of the analytical logic,
- documentation of software versions and the computational environment, and
- a record of analytical and curation decisions (thresholds and settings in `config/`, KD group assignments in `metadata/`).
Inputs are described in the manuscript Methods. The scripts are not packaged for one-click execution.

---
 
## Notes on manual curation
 
Some inputs and tables were manually harmonized prior to analysis or inclusion. These operations included:
 
- reclassifying PHASTER prophage labels against vOTU references
- manually classifying BF (BLAST family) families from TAfinder output
- reviewing the outputs of S17A, S17B and S17C before the next step (removal of hypothetical and essential genes, gene categorization, paralog review)
- standardizing categorical labels and identifiers
- reformatting and harmonizing matrix structures for consistency
- renaming column headers
These curation steps are documented in the manuscript and its Supplementary Information. An automated re-run of S17A–S17D reproduces the pipeline but not the manual review between steps.
