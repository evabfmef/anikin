# anikin
Code and analysis scripts for **Kin discrimination in *Bacillus subtilis* is predicted by average nucleotide identity but shaped by mobile genetic elements**; Stare et al., 2026 (to be updated, in prep.)

---
## Overview
This repository contains analysis code and documentation supporting the manuscript:

Title: Kin discrimination in *Bacillus subtilis* is predicted by average nucleotide identity but shaped by mobile genetic elements

Authors: to be updated

Journal: to be updated

DOI: to be updated

This study examines the genomic basis of kin discrimination in *Bacillus subtilis* by combining comparative genomics, swarming assays (850 pairwise combinations), and machine learning. We identify a species-wide ANI gap at ~99.5% that delineates socially compatible from incompatible populations, and show that genome-wide relatedness alone does not predict social compatibility among closely related strains. Instead, accessory determinants act as the additional drivers of discrimination outcomes.

---

## Abstract

### Background
Social interactions in bacteria can generate cooperative subpopulations, yet the genomic principles defining social compatibility within natural populations remain poorly understood. In *Bacillus subtilis*, kin discrimination (KD) during swarming enables strains to distinguish related from unrelated partners, but its genomic basis remains unresolved. Using 313 high-quality genomes, we identified several genomic
determinants of KD through group-specific feature association and supervised machine-learning classification.

### Results
Average nucleotide identity (ANI; computed with fastANI), core-genome phylogeny, and accessory gene content converged on comparable strain groupings. ANI analysis revealed a pronounced underrepresentation of genome pairs at 99.0–99.5% ANI, both locally and species-wide — an intraspecific "ANI gap." Phenotypically, ~99.5% ANI marked the kin boundary, yet this threshold proved necessary but not sufficient: 26.3%
of strain pairs above it showed incompatibility, including lysis between strains sharing >99.9% ANI. Three analytical frameworks converged on a multi-layered KD architecture — recognition through cell-surface gene cassettes, antagonism via LXG toxin diversification and SPβ prophage mosaicism, and communication/exclusion through variation in the conjugative element ICEBs1, quorum sensing, and spore-coat composition.

### Conclusion
The intraspecific ANI gap coincides with the functional boundary of social compatibility in *B. subtilis*, and incompatibility between highly related strains is driven less by core-genome divergence than by asymmetric carriage of mobile genetic elements. ANI is established as a practical predictor of kin group boundaries, the ANI gap is confirmed as a species-wide phenomenon, and accessory-genome mosaicism —
particularly SPβ prophage variation, LXG polymorphic toxin/antitoxin variants, and asymmetric carriage of mobile-element cargo (PBSX, ICEBs1, skin element) — is identified as the primary determinant overriding genome-wide relatedness.

## Key Findings
- Kin discrimination groups in *Bacillus subtilis* closely mirror genome-wide population structure and are governed by a multi-layered, modular accessory-genome architecture rather than by single conserved determinants.
- A conserved species-wide "ANI gap" at ~99.5% ANI separates socially compatible and incompatible *B. subtilis* populations, and fastANI reliably predicts kin group association.
- Closely related strains can remain socially incompatible despite sharing >99.9% genome identity, with incompatibility driven primarily by asymmetric carriage of SPβ prophages, LXG toxin systems, and other mobile genetic elements.

---
 
## Repository structure
 
```
anikin/
├── README.md                     # this file
├── scripts/
│   ├── README.md                 # script index and notes
│   └── S1–S18 ...                # analysis and figure-generation scripts
├── docs/
│   ├── data_overview.md          # datasets and curation notes
│   └── software_versions.md      # tool, database, and environment versions
└── environment/
    ├── README.md                 # environment notes
    ├── R_packages.txt            # R packages used by the R scripts
    └── requirements.txt          # Python packages used by the Python scripts
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
