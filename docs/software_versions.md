# Software versions

## Purpose

This document records the main software tools, databases, and computational environments used in the study *"Kin discrimination in Bacillus subtilis is predicted by average nucleotide identity but shaped by mobile genetic elements"* (Stare et al.).

Versions are those used in this study. Tools without a version below were not version-pinned in the original workflow; where known, the version is given. Python and R package versions are listed separately in `environment/requirements.txt` and `environment/R_packages.txt`.

---

## Sequencing data processing

**Dorado** (Nanopore basecalling; dna_r10.4.1_e8.2_400bps_hac@v4.2.0 model)
- Version: v0.4.1

**Trimmomatic** (Illumina read processing)
- Version: v38.1 / v0.39

**NanoFilt** (Nanopore sequencing data processing)
- Version: v2.8.0

**Canu** (PacBio and Nanopore read processing)
- Version: v2.1.1 / v2.2

**NanoPack** (package of tools for Nanopore sequencing data evaluation)
- NanoStat
  - Version: v1.6.0
- NanoPlot
  - Version: v1.40.0
- NanoQC
  - Version: v0.9.4

**NCBI Datasets CLI** (retrieval of public genomes)
- Version: v18.15.0

---

## Genome assembly

**Unicycler** (Galaxy Europe platform)
- Version: v0.5.0

**QUAST** (genome assembly evaluation)
- Version: v5.0.2

---

## Genome quality assessment

**BUSCO** (completeness; bacillales_odb10 lineage dataset)
- Version: v6.0.0

**CheckM** (completeness and contamination)
- Version: v1.2.4

---

## Functional annotation

**Prokka** (genome annotation, custom reference database)
- Version: v1.14.6

**antiSMASH** (biosynthetic gene cluster prediction)
- Version: v6.1.1

**PHASTER** (prophage prediction)
- Version: no version number; accessed via PHASTER server in 2024

**TAfinder** (toxin–antitoxin system identification)
- Version: v2.0 (database: TADB v3.0)

**PredicTF** (transcription factor prediction)
- Version: no version pinned; accessed in October 2025

---

## Average nucleotide identity and pangenome

**fastANI** (average nucleotide identity)
- Version: v1.33

**Roary** (pangenome analysis; Galaxy Europe platform)
- Version: v3.13.0

---

## Phylogenetic analysis

**RAxML** (maximum-likelihood core-genome phylogeny; Galaxy Europe platform)
- Version: v8.2.12

**MEGA11** (building phylogenetic trees and dendrograms)
- Version: v11.0.13

---

## Sequence comparison and curation

**BLAST+** (sequence similarity searches; dc-megablast, blastp, makeblastdb)
- Version: 2.12.0+

**MUSCLE** (multiple sequence alignment of the yxi/deaD-yxxF region; run within Geneious, PPP algorithm)
- Version: v5.1

**Geneious Prime** (sequence extraction and manual inspection)
- Version: 2025.2.2

---

## Statistical analysis and scripting

**R**
- Version: 4.5.0
- Key packages: see `environment/R_packages.txt`

**Python**
- Version: 3.10 / 3.12
- Key libraries: see `environment/requirements.txt`

---

## Databases and reference resources

- **SubtiWiki** v5 (gene annotation; Elfmann et al., 2025)
- **TADB** v3.0 (toxin-antitoxin reference, used by TAfinder)
- **NCBI RefSeq** (reference proteome and annotation)
- **Prophage vOTU references** in Štefanič et al. (2025), Zenodo 13834859
