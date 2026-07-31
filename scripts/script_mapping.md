\# Script I/O map



Inputs, outputs, and dependencies for the analysis scripts (S1–S18). Script filenames use canonical names (see `README.md` for the script index).



| script | input | output | depends\_on | Notes |

| --- | --- | --- | --- | --- |

| S1\_Prokka\_annotation.sh | Genome assemblies (\*.fasta) + custom Prokka protein DB (manual\_database.fasta); requires Prokka | Per-strain Prokka annotations (GFF, GBK, FAA) in annot\_output/<strain>/ | - |  |

| S2\_fastANI\_analysis.R | fastANI classical output (TSV, no header: query, reference, ani, coverage, fragments) | ANI heatmaps, distribution overlay, descriptive stats, gap-test stats, clustered matrix (.xlsx) | - | (line 735) "highres=TRUE" is recommended for 313 strain analysis |

| S3\_ANI\_density\_overlay.R | Two fastANI classical TSVs (40-strain + 313-strain) | Kernel-density overlay plot (SVG) in figures/ | - |  |

| S4\_pangenome\_plots.R | Roary .Rtab files (new/conserved/pan/unique genes + blast identity frequencies) | Pangenome accumulation + BLAST-identity plots | - |  |

| S5A\_conservative\_gene\_clustering.py | mega\_matrix.xlsx, locus\_tag\_mapping.xlsx, strain\_list.xlsx, GBK folder | Per-gene JSON + 01\_conservative\_clustering\_summary.xlsx | - |  |

| S5B\_assembly\_artifact\_detection.py | Per-gene JSON files from S5A | 02\_assembly\_artifact\_detection\_summary.xlsx | S5A |  |

| S5C\_merge\_variants.py | mega\_matrix.xlsx + S5A summary + S5B summary | mega\_matrix\_merged.xlsx (+ merge mapping/report) | S5A, S5B |  |

| S6\_PHASTER\_matrix\_generation.R | phaster\_results.xlsx (sheet: All Phaster Results) | Prophage presence/absence matrix (.xlsx) | - |  |

| S7\_phaster\_to\_votu\_mapping.py | PHASTER output dir + vOTU reference FASTAs; requires BLAST+ | Prophage->vOTU mapping TSV, crosstab TSV, strain x vOTU matrix TSV | - |  |

| S8\_TA\_systems\_matrix\_generation.R | TA\_pairs\_results.xlsx (sheet: Split\_product) | TA family presence/absence matrix (.xlsx) | - |  |

| S9\_BGC\_matrix\_generation.R | antismash\_results.xlsx (sheets: Sheet 1 + chromosome\_length) | BGC presence/absence matrix (.xlsx) | - |  |

| S10A\_GEM\_analysis\_part1.py | Reference model iBB1018.xml, per-strain GBK, StrainInformation.xlsx; requires BLAST+ | BBH results, ortho\_matrix.csv, geneIDs\_matrix.csv | - |  |

| S10B\_GEM\_analysis\_part2.py | iBB1018.xml + S10A ortho/geneID matrices | Per-strain models GEM/Models/<strain>.json | S10A |  |

| S10C\_GEM\_analysis\_part3.py | S10B \*.json models | reaction\_matrix.xlsx (binary reaction presence/absence) | S10B |  |

| S11A\_batch\_predicTF.sh | Protein \*.fa/\*.faa files; requires PredicTF | Per-genome file.out.mapping.potential.TF | - |  |

| S11B\_TF\_matrix\_generation.R | Folder of PredicTF \*.potential.TF files (from S11A) | TF\_matrix.xlsx | S11A |  |

| S12\_gene\_content\_clustering.R | Mega\_matrix.xlsx (sheet: Filtered\_Genes, binary matrix) | UPGMA tree gene\_content\_UPGMA\_rescaled.nwk | - |  |

| S13A\_strain\_compatibility\_analysis.R | strain\_interactions.xlsx (cols: S1, S2, Phenotype\_interaction, Group\_S1, Group\_S2, fastANI); paths/thresholds/phenotype defs read from config\_S13.yml | ANI threshold summary stats + fastANI distribution plot (Fig 2C); nonconforming pairs stratified into strong/moderate/borderline classes | - |  |

| S13B\_threshold\_validation.R | strain\_interactions.xlsx (same file as S13A); config\_S13.yml | Threshold uncertainty + predictive validation: bootstrap CIs for per-class minimum ANI (strain-level resampling), sensitivity sweep 99.3-99.8%, ROC (AUC + DeLong CI) \& precision-recall, leave-group-out validation (threshold\_bootstrap\_ci.csv, threshold\_sensitivity\_sweep.csv, threshold\_roc\_pr.csv, threshold\_leave\_group\_out.csv) | S13A |  |

| S14\_KD\_group\_validation.py | Strain\_interactions.xlsx (cols: S1, S2, Phenotype\_interaction, Group\_S1, Group\_S2) | KD\_Validation\_with\_Subgroups.xlsx | - |  |

| S15A\_group\_specific\_features.py | megamatrix\_filtered\_strains.xlsx (binary matrix; group membership hard-coded) | features\_analysis\_results.xlsx | - |  |

| S15B\_KD\_variant\_validation.py | Megamatrix .xlsx + Features\_analysis\_67\_output.xlsx + GBK folder; requires Biopython | variant\_validation\_results.xlsx + corrected\_features\_analysis.xlsx | - |  |

| S15C\_paralog\_check.py | Features\_analysis\_67\_output.xlsx + GBK folder (67\_strains\_gbks/) | Paralog-check results (.xlsx) | - |  |

| S15D\_paralog\_annotation.py | Manually reviewed: 67\_corrected\_features\_analysis\_reviewed.xlsx, 67\_paralog\_check\_results\_reviewed.xlsx + 67\_mega\_matrix\_full.xlsx | Paralog-annotated features table | S15C (manual review) |  |

| S15E\_generate\_sankey\_3level\_weighted.py | 39\_KD\_features\_summary\_9kat\_p\_PRESENT.csv (semicolon-separated; cols: Kin group, Gene name, Functional category, Category) | 3-level weighted Sankey diagram | - |  |

| S16\_ML\_classification.py | 67\_megamatrix\_for\_ML.xlsx (or 39-strain version; binary matrix) | Trained models \*.joblib + feature\_patterns\_by\_kin\_group\_NB.xlsx, feature\_importances.xlsx, robust\_candidates.xlsx | - |  |

| S17\_feature\_enrichment\_and\_strain\_centric\_analysis.R | strain\_interactions.xlsx + Filtered\_presence\_absence\_78strains.xlsx | Enrichment tables + figures | - |  |

| S18\_generate\_chromosome\_map.py | strain\_centric\_recurrent\_candidates\_annotated\_reviewed.xlsx (manually reviewed) | Chromosome map SVG | manual review |  |

