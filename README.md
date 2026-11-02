# anikin

Code and analysis scripts for **Genomic boundaries mark social boundaries in *Bacillus subtilis***; Stare et al., 2026 (in preparation).

---

## Overview

This repository contains analysis code and documentation supporting the manuscript:

- **Title:** Genomic boundaries mark social boundaries in *Bacillus subtilis*
- **Authors:** Eva Stare, Manca Mršol, Robson Parmezan Bonidia, Anna Dragoš, Jasna Kovac, Ákos T. Kovács, Ines Mandic-Mulec, Ulisses Nunes da Rocha, Polonca Stefanic
- **Journal:** The ISME Journal
- **DOI:** to be updated

This study examines how genomic relatedness structures social compatibility in natural populations of *Bacillus subtilis*, combining comparative genomics of 313 genomes with 850 pairwise swarming interactions, enrichment analysis and machine learning. A ~99.5% average nucleotide identity (ANI) boundary separates strains that discriminate against one another and coincides with an intraspecific ANI gap, while discrimination among highly related strains is associated with modular variation in the accessory genome.

The abstract and key findings will be added at publication.

---

## Repository structure

```
anikin/
├── README.md                     # this file
├── scripts/
│   ├── README.md                 # script index and notes
│   └── S1–S20 ...                # analysis and figure-generation scripts
├── config/
│   ├── config_S15.yml            # settings for S15A, S15B
│   ├── config_S16_S19.yml        # settings for S17A–S17D, S19
│   └── kd_config.py              # loads config_S16_S19.yml and metadata/ for S17A–S17D
├── metadata/
│   ├── accessions.tsv            # 313 genomes: NCBI accessions and genome file names
│   ├── kd_groups_39.tsv          # KD group membership, 39-strain panel
│   ├── kd_groups_67.tsv          # KD group membership, 67-strain panel
│   ├── kd_groups_representatives_39.tsv
│   └── kd_groups_representatives_67.tsv
├── docs/
│   ├── data_overview.md          # datasets, repository contents and curation notes
│   └── software_versions.md      # tool, database and environment versions
└── environment/
    ├── environment.yml           # conda environment (all Python scripts except S11B, S11C)
    ├── requirements.txt          # same packages for pip
    ├── requirements_GEM.txt      # separate environment for S11B, S11C
    └── R_packages.txt            # R packages used by the R scripts
```

---

## Reproducibility scope

This repository provides the analysis scripts and documents the software versions and environment underlying the manuscript. It is intended to make the analysis steps transparent and inspectable. Inputs are described in the manuscript Methods; the scripts are not packaged for one-click execution. Scripts are run from the repository root (see `scripts/README.md`).

#### Included in this repository

- Analysis and figure-generation scripts (S1–S20)
- Analysis settings (`config/`) and strain-level tables: genome accessions and KD group assignments (`metadata/`)
- Package lists (`environment/`)
- Software, database and environment versions (`docs/software_versions.md`)
- Datasets, provenance and curation notes (`docs/data_overview.md`)

The main and supplementary figures and the final curated supplementary tables are published with the article and are not duplicated here.

#### Not included

- Raw sequencing reads
- Genome assemblies and annotations
- Raw outputs of bioinformatics tools and intermediate files
- Manuscript figures, supplementary figures and supplementary tables (available with the article)

---

## Requirements

The workflow uses bioinformatics tools (e.g. fastANI, Roary, Prokka, antiSMASH) together with Python (e.g. scikit-learn, pandas, Biopython) and R (e.g. tidyverse, ggplot2, ape). Python scripts run in the environment in `environment/environment.yml`, except S11B and S11C, which need the separate environment in `environment/requirements_GEM.txt`. R packages are listed in `environment/R_packages.txt`; the complete tool, database and version list is in `docs/software_versions.md`.

---

## Citation

To be updated at publication.