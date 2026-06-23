# Scripts
 
## Overview
 
This directory contains the analysis and figure-generation scripts (S1–S18) underlying the study *"Average nucleotide identity predicts kin discrimination, but mobile genetic elements determine its outcome in Bacillus subtilis"* (Stare et al.). They implement the analytical steps described in the Methods, from genome annotation and average nucleotide identity (ANI), pangenome analysis, feature-matrix construction, machine-learning classification, enrichment study, and figure generation.
 
Scripts are numbered in approximate order of use. Multi-part steps share a number with letter suffixes (e.g. S5A, S5B, S5C).
 
## Notes
 
- Input and output paths are defined at the top of each script and may need to be adjusted to a local setup.
- Inputs are described in the manuscript Methods; primary sequence data and intermediate files are not redistributed here (see `docs/data_overview.md`).
- Software and package versions are recorded in `docs/software_versions.md` and `environment/`.

## Reproducibility statement
 
These scripts document and implement the analytical strategy described in the Methods section of the manuscript. They are provided for transparency and inspection of the analytical logic, rather than as a packaged, reproducible pipeline.
