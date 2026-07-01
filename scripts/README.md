# Scripts
 
## Overview
 
This directory contains the analysis and figure-generation scripts (S1–S18) underlying the study *"Kin discrimination in Bacillus subtilis is predicted by average nucleotide identity but shaped by mobile genetic elements"* (Stare et al.). They implement the analytical steps described in the Methods, from genome annotation and average nucleotide identity (ANI), pangenome analysis, feature-matrix construction, machine-learning classification, enrichment study, and figure generation.
 
Scripts are numbered in approximate order of use. Multi-part steps share a number with letter suffixes (e.g. S5A, S5B, S5C).

## Script index
 
| Script | Language | Description |
| --- | --- | --- |
| S1 | bash | Genome annotation with Prokka using a custom reference database |
| S2 | R | fastANI analysis |
| S3 | R | ANI density overlay plot across strain panels |
| S4 | R | Roary pangenome summary plots |
| S5A | Python | Conservative gene-variant clustering |
| S5B | Python | Assembly-artifact detection |
| S5C | Python | Merge gene variants into a reduced presence/absence matrix |
| S6 | R | PHASTER prophage presence/absence matrix generation |
| S7 | Python | PHASTER-to-vOTU reclassification via BLAST against vOTU references |
| S8 | R | Toxin–antitoxin (TAfinder) systems matrix generation |
| S9 | R | Biosynthetic gene cluster (antiSMASH) matrix generation |
| S10A | Python | GEM analysis: bidirectional BLASTp and homology matrix |
| S10B | Python | GEM analysis: draft strain-specific metabolic models |
| S10C | Python | GEM analysis: reaction presence/absence matrix |
| S11A | bash | Batch PredicTF transcription factor prediction |
| S11B | R | Transcription factor presence/absence matrix generation |
| S12 | R | Accessory gene-content clustering |
| S13 | R | Strain compatibility analysis |
| S14 | Python | KD group validation and subgroup generation |
| S15A | Python | Identification of KD group-specific genomic features |
| S15B | Python | Validation of KD group-specific variants |
| S15C | Python | Paralog check for candidate features |
| S15D | Python | Paralog annotation |
| S15E | Python | Three-level weighted Sankey diagram of group-specific features |
| S16 | Python | Machine-learning classification of kin discrimination group membership |
| S17 | R | Feature enrichment and strain-centric case–control analysis |
| S18 | Python | Chromosome map of recurrent candidate genes |
 

## Notes
 
- Input and output paths are defined at the top of each script and may need to be adjusted to a local setup.
- Inputs are described in the manuscript Methods; primary sequence data and intermediate files are not redistributed here (see `docs/data_overview.md`).
- Software and package versions are recorded in `docs/software_versions.md` and `environment/`.

## Reproducibility statement
 
These scripts document and implement the analytical strategy described in the Methods section of the manuscript. They are provided for transparency and inspection of the analytical logic, rather than as a packaged, reproducible pipeline.
