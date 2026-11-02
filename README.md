# anikin
Code and analysis scripts for **Genomic boundaries mark social boundaries in Bacillus subtilis*; Stare et al., 2026 (to be updated, in prep.)

---
## Overview
This repository contains analysis code and documentation supporting the manuscript:
 
Title: Genomic boundaries mark social boundaries in *Bacillus subtilis*
 
Authors: Eva Stare, Manca Mršol, Robson Parmezan Bonidia, Anna Dragoš, Jasna Kovac, Ákos T. Kovács, Ines Mandic-Mulec, Ulisses Nunes da Rocha, Polonca Stefanic
 
Journal: The ISME Journal: Multidisciplinary Journal of Microbial Ecology
 
DOI: to be updated
 
This study examines the genomic basis of kin discrimination in *Bacillus subtilis* by combining comparative genomics, swarming assays (850 pairwise combinations), and machine learning. We identify a phenotypic ANI boundary at ~99.5% that delineates socially compatible from incompatible populations, and show that genome-wide relatedness alone does not predict social compatibility among closely related strains. Instead, accessory determinants act as the additional drivers of discrimination outcomes.
 
---
 
## Abstract
To be updated at publication
 
### Results
To be updated at publication
 
### Conclusion
To be updated at publication
 
## Key Findings
To be updated at publication

---
 
## Repository structure
 
```
anikin/
├── README.md                     # this file
├── scripts/
│   ├── README.md                 # script index and notes
│   ├── script_mapping.md         # per-script input/output/dependency map
│   └── S1–S18 ...                # analysis and figure-generation scripts
├── docs/
│   ├── data_overview.md          # datasets and curation notes
│   └── software_versions.md      # tool, database, and environment versions
└── environment/
    ├── README.md                 # environment notes
    ├── environment.yml           # conda environment
    ├── R_packages.txt            # R packages used by the R scripts
    ├── requirements.txt          # Python packages used by the Python scripts
    └── requirements_GEM.txt      # Python packages for the GEM analysis (S10)
```
 
---
  
## Reproducibility Scope
This repository provides the analysis scripts and documents the software versions and environment underlying the manuscript. It is intended to make the analysis steps transparent and inspectable. Inputs are described in the manuscript Methods section; the scripts are not packaged for one-click execution.
 
#### Included in this repository
- Analysis and figure-generation scripts (S1–S18)
- Package lists (`environment/`)
- Software, database, and environment versions (`docs/software_versions.md`)
- Datasets, provenance, and curation notes (`docs/data_overview.md`)
Note: The main and supplementary figures and final curated supplementary tables are published with the article and are not duplicated here; see the Supplementary Information accompanying the manuscript.
 
#### Not included
- Raw sequencing reads
- Genome assembly workflows
- Genome assemblies and any intermediate files
- Manuscript figures, supplementary figures and supplementary tables (available with the article)
---
 
## Requirements
 
In brief, the workflow uses bioinformatics tools (e.g. fastANI, Roary, Prokka, antiSMASH) together with Python (e.g. scikit-learn, pandas, BioPython) and R (e.g. tidyverse, ggplot2, ape). Full package lists are in `environment/`; the complete tool, database, and version list is in `docs/software_versions.md` and the manuscript Methods.
 
---
## Citation
to be updated