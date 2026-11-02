# Scripts

## Overview

This directory contains the analysis and figure-generation scripts (S1–S20) underlying the study *"Genomic boundaries mark social boundaries in Bacillus subtilis"* (Stare et al.). They implement the analytical steps described in the Methods, from genome annotation and average nucleotide identity (ANI), pangenome and phylogenetic concordance analysis, feature-matrix construction, ANI–kin discrimination (KD) statistics, machine-learning classification and enrichment analysis to figure generation.

Scripts are numbered in approximate order of use. Multi-part steps share a number with letter suffixes (e.g. S5A, S5B, S5C).

## Script index

### Genome annotation, ANI and pangenome

| Script | Language | Description |
| --- | --- | --- |
| S1 | bash | Genome annotation with Prokka using a custom reference database |
| S2 | R | fastANI analysis |
| S3 | R | ANI density overlay plot across strain panels |
| S4 | R | Roary pangenome summary plots |
| S5A | Python | Conservative gene-variant clustering |
| S5B | Python | Assembly-artifact detection |
| S5C | Python | Merge gene variants into a reduced presence/absence matrix |
| S6 | Python | Concordance of ANI genomovars (≥99.5% ANI) with the core-genome phylogeny |

### Feature-matrix construction

| Script | Language | Description |
| --- | --- | --- |
| S7 | R | PHASTER prophage presence/absence matrix generation |
| S8 | Python | PHASTER-to-vOTU reclassification via BLAST against vOTU references |
| S9 | R | Toxin–antitoxin (TAfinder) systems matrix generation |
| S10 | R | Biosynthetic gene cluster (antiSMASH) matrix generation |
| S11A | Python | GEM analysis: bidirectional BLASTp and homology matrix |
| S11B | Python | GEM analysis: draft strain-specific metabolic models |
| S11C | Python | GEM analysis: reaction presence/absence matrix |
| S12A | bash | Batch PredicTF transcription factor prediction |
| S12B | R | Transcription factor presence/absence matrix generation |
| S13 | R | Accessory gene-content clustering |

### ANI and kin discrimination

| Script | Language | Description |
| --- | --- | --- |
| S14 | Python | KD phenotype heatmap (Fig. 3B) and ANI–KD statistics: correlation, strain-label permutation test and Bayesian ordinal mixed-effects model |
| S15A | R | Strain compatibility analysis |
| S15B | R | fastANI kin-threshold validation: strain-level bootstrap, sensitivity sweep, ROC/PR and leave-group-out validation |
| S16 | Python | KD group validation and subgroup generation |

### KD group-specific features, classification and enrichment

| Script | Language | Description |
| --- | --- | --- |
| S17A | Python | Identification of KD group-specific genomic features |
| S17B | Python | Validation of KD group-specific variants |
| S17C | Python | Paralog check for candidate features |
| S17D | Python | Paralog annotation |
| S17E | Python | Three-level weighted Sankey diagram of group-specific features |
| S18 | Python | Machine-learning classification of KD group membership |
| S19 | R | Feature enrichment and strain-centric case–control analysis |
| S20 | Python | Chromosome map of recurrent candidate genes |

## Notes

- Run scripts from the repository root: configuration, metadata and input paths are relative to it.
- Input and output paths are defined at the top of each script, in `config/` or passed as command-line arguments, and may need to be adjusted to a local setup.
- S15A and S15B share `config/config_S15.yml` (paths, phenotype scale, threshold-validation settings), so the two scripts use identical settings. From the repository root, pass `--config=config/config_S15.yml`.
- S17A–S17D and S19 read their settings from `config/config_S16_S19.yml`. S17A–S17D load it, together with the KD group and accession tables in `metadata/`, through the helper module `config/kd_config.py`; select the strain panel with `--analysis=39` or `--analysis=67`. For S19, pass `--config=config/config_S16_S19.yml`.
- S6 and S14 read directly from the Supplementary Tables workbook (S6: sheet 8; S14: sheets 5, 9 and 10).
- S11B and S11C run in a separate Python environment (`environment/requirements_GEM.txt`); all other Python scripts use `environment/environment.yml`.
- Inputs are described in the manuscript Methods; primary sequence data and intermediate files are not redistributed here (see `docs/data_overview.md`).
- Software and package versions are recorded in `docs/software_versions.md` and `environment/`.

## Reproducibility statement

These scripts document and implement the analytical strategy described in the Methods section of the manuscript. They are provided for transparency and inspection of the analytical logic, rather than as a packaged, reproducible pipeline.
