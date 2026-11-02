# Script I/O map

Inputs, outputs and dependencies for the analysis scripts (S1–S20). Script filenames use canonical names (see `README.md` for the script index).

| Script | Input | Output | Depends on | Notes |
| --- | --- | --- | --- | --- |
| S1_Prokka_annotation.sh | Genome assemblies (*.fasta) + custom Prokka protein DB (manual_database.fasta); requires Prokka | Per-strain Prokka annotations (GFF, GBK, FAA) in annot_output/&lt;strain&gt;/ | – | |
| S2_fastANI_analysis.R | fastANI classical output (TSV, no header: query, reference, ani, coverage, fragments) | ANI heatmaps, distribution overlay, descriptive stats, gap-test stats, clustered matrix (.xlsx) | – | Set `highres = TRUE` in the main entry-point call (section 9) for the 313-strain analysis |
| S3_ANI_density_overlay.R | Two fastANI classical TSVs (40-strain + 313-strain) | Kernel-density overlay plot (SVG) in figures/ | – | |
| S4_pangenome_plots.R | Roary .Rtab files (new/conserved/pan/unique genes + BLAST identity frequencies) | Pangenome accumulation + BLAST-identity plots | – | |
| S5A_conservative_gene_clustering.py | mega_matrix.xlsx, locus_tag_mapping.xlsx, strain_list.xlsx, GBK folder | Per-gene JSON + 01_conservative_clustering_summary.xlsx | – | |
| S5B_assembly_artifact_detection.py | Per-gene JSON files from S5A | 02_assembly_artifact_detection_summary.xlsx | S5A | |
| S5C_merge_variants.py | mega_matrix.xlsx + S5A summary + S5B summary | mega_matrix_merged.xlsx (+ merge mapping/report) | S5A, S5B | |
| S6_genomovar_concordance.py | Core-genome Newick tree with bootstrap support (core313.nwk) + Supplementary Tables workbook (sheet 8: genomovars); optional aliases JSON; requires Python ≥ 3.10, openpyxl | genomovar_results.csv, strain_mapping.csv, summary.json (concordance %, software versions, SHA-256 checksums) in output/S6/ | – | Supplementary Methods 4 |
| S7_PHASTER_matrix_generation.R | phaster_results.xlsx (sheet: All Phaster Results) | Prophage presence/absence matrix (.xlsx) | – | |
| S8_phaster_to_votu_mapping.py | PHASTER output dir + vOTU reference FASTAs; requires BLAST+ | Prophage→vOTU mapping TSV, crosstab TSV, strain × vOTU matrix TSV | – | |
| S9_TA_systems_matrix_generation.R | TA_pairs_results.xlsx (sheet: Split_product) | TA family presence/absence matrix (.xlsx) | – | |
| S10_BGC_matrix_generation.R | antismash_results.xlsx (sheets: Sheet 1 + chromosome_length) | BGC presence/absence matrix (.xlsx) | – | |
| S11A_GEM_analysis_part1.py | Reference model iBB1018.xml, per-strain GBK, StrainInformation.xlsx; requires BLAST+ | BBH results, ortho_matrix.csv, geneIDs_matrix.csv | – | |
| S11B_GEM_analysis_part2.py | iBB1018.xml + S11A ortho/geneID matrices | Per-strain models GEM/Models/&lt;strain&gt;.json | S11A | |
| S11C_GEM_analysis_part3.py | S11B *.json models | reaction_matrix.xlsx (binary reaction presence/absence) | S11B | |
| S12A_batch_predicTF.sh | Protein *.fa/*.faa files; requires PredicTF | Per-genome file.out.mapping.potential.TF | – | |
| S12B_TF_matrix_generation.R | Folder of PredicTF *.potential.TF files (from S12A) | TF_matrix.xlsx | S12A | |
| S13_gene_content_clustering.R | Mega_matrix.xlsx (sheet: Filtered_Genes, binary matrix) | UPGMA tree gene_content_UPGMA_rescaled.nwk | – | |
| S14_KD_ANI_heatmap_and_statistics.py | Supplementary Tables workbook (sheet 5: ANI matrix, 313 strains; sheet 9: candidate KD groups; sheet 10: 850 tested pairs); requires pandas, numpy, scipy, matplotlib, openpyxl, pymc, arviz | KD heatmap (PNG/PDF/SVG), permutation null distribution CSV, ordinal model summary CSV, posterior predictive frequencies CSV, model trace (.nc) in output/S14/ | – | Fig. 3B; Supplementary Methods 8. `--skip-model` runs heatmap and permutation test only. Published analysis run in Google Colab |
| S15A_strain_compatibility_analysis.R | strain_interactions.xlsx (cols: S1, S2, Phenotype_interaction, Group_S1, Group_S2, fastANI); paths/thresholds/phenotype definitions read from config_S15.yml | Within-group summary, fastANI threshold summary, Pearson correlation, nonconforming pairs (strong/moderate/borderline CSVs), fastANI distribution plot (Fig. 2C), summary.txt | – | |
| S15B_threshold_validation.R | strain_interactions.xlsx (same file as S15A); config_S15.yml | Threshold uncertainty + predictive validation: bootstrap CIs for per-class minimum ANI (strain-level resampling), sensitivity sweep 99.3–99.8%, ROC (AUC + DeLong CI) and precision-recall, leave-group-out validation (threshold_bootstrap_ci.csv, threshold_sensitivity_sweep.csv, threshold_roc_pr.csv, threshold_leave_group_out.csv, threshold_validation.png, threshold_validation_summary.txt) | S15A | |
| S16_KD_group_validation.py | Strain_interactions.xlsx (cols: S1, S2, Phenotype_interaction, Group_S1, Group_S2) | KD_Validation_with_Subgroups.xlsx | – | |
| S17A_group_specific_features.py | megamatrix_filtered_strains.xlsx (binary matrix; group membership hard-coded) | features_analysis_results.xlsx | – | |
| S17B_KD_variant_validation.py | Megamatrix .xlsx + Features_analysis_67_output.xlsx + GBK folder; requires Biopython | variant_validation_results.xlsx + corrected_features_analysis.xlsx | – | |
| S17C_paralog_check.py | Features_analysis_67_output.xlsx + GBK folder (67_strains_gbks/) | Paralog-check results (.xlsx) | – | |
| S17D_paralog_annotation.py | Manually reviewed: 67_corrected_features_analysis_reviewed.xlsx, 67_paralog_check_results_reviewed.xlsx + 67_mega_matrix_full.xlsx | Paralog-annotated features table | S17C (manual review) | |
| S17E_generate_sankey_3level_weighted.py | 39_KD_features_summary_9kat_p_PRESENT.csv (semicolon-separated; cols: Kin group, Gene name, Functional category, Category) | Three-level weighted Sankey diagram | – | |
| S18_ML_classification.py | 67_megamatrix_for_ML.xlsx (or 39-strain version; binary matrix) | Trained models *.joblib + feature_patterns_by_kin_group_NB.xlsx, feature_importances.xlsx, robust_candidates.xlsx | – | |
| S19_feature_enrichment_and_strain_centric_analysis.R | strain_interactions.xlsx + Filtered_presence_absence_78strains.xlsx | Enrichment tables + figures | – | |
| S20_generate_chromosome_map.py | strain_centric_recurrent_candidates_annotated_reviewed.xlsx (manually reviewed) | Chromosome map SVG | Manual review | |