# =============================================================================
# config_S16_S19.yml
#
# Configuration for S17A-S17D and S19. (S16 is self-contained; its inputs and
# the pinned minimum vertex cover are defined at the top of the script.)
#
# All thresholds, reference-strain choices and inferential cut-offs live here
# so that biological decisions are reviewable without reading Python or R
# internals. Strain-level data (KD group memberships, representative strains,
# genome accessions) lives in TSV files under metadata/, listed in the
# `metadata:` block below.
#
# Companion file: config_S15.yml (S15A / S15B).
# =============================================================================

# --- Analysis panel -----------------------------------------------------------
# Selects which strain panel the S17 series operates on. Overridable on the
# command line with --analysis=39.
analysis: "67"

# --- Metadata tables ----------------------------------------------------------
# Strain-level biological data. {analysis} is substituted at run time.
metadata:
  kd_groups: "metadata/kd_groups_{analysis}.tsv"
  kd_groups_representatives: "metadata/kd_groups_representatives_{analysis}.tsv"
  # Strain name -> genome file name (public genomes are stored under their
  # assembly accession). Also lists the NCBI accession of every genome.
  accessions: "metadata/accessions.tsv"

# --- Phenotype coding ---------------------------------------------------------
# Shared with config_S15.yml; kept here so S19 is self-contained.
#   1.0  = strong kin (merge)        0.25 = weak non-kin (discrimination)
#   0.5  = weak kin (merge)          0.0  = strong non-kin (discrimination)
phenotype:
  kin_values: [1.0, 0.5]
  discrimination_values: [0.25, 0.0]
  kin_cutoff: 0.5
  # Free-text phenotype synonyms accepted by S19's parser, mapped to values.
  text_aliases:
    Merging: 1.0
    Contact: 0.5
    Thick: 0.25
    Lysis: 0.0
  # Display labels for phenotype categories (S19).
  category_labels:
    1.0: "Merging"
    0.5: "Contact"
    0.25: "Thick boundary"
    0.0: "Lysis"

# =============================================================================
# S17A - group-specific features
# =============================================================================
s17a:
  input_file_67: "data/kd/inputs/67_matrix_strains_as_rows.tsv"
  input_file_39: "data/matrices_40PS/40PS_matrix_genes.tsv"
  output_dir: "data/kd/output/{analysis}"
  output_file: "features_analysis_results.xlsx"
  summary_file: "group_size_summary.xlsx"

  # Group size classes. The exact-pattern criterion is descriptive, and its
  # interpretation depends on group size: in a singleton it reduces to private
  # accessory content, and in a pair it is satisfied readily by chance.
  size_classes:
    singleton_max: 1
    pair_max: 2
    # anything larger is MULTI

# =============================================================================
# S17B - variant validation by pairwise protein alignment
# =============================================================================
s17b:
  base_dir: "data/kd/"
  gbk_dir: "data/raw/gbk"
  megamatrix_67: "data/kd/inputs/67_matrix_strains_as_rows.tsv"
  megamatrix_39: "data/kd/inputs/39_matrix_strains_as_rows.tsv"
  features_file: "features_analysis_results.xlsx"   # S17A output
  output_dir: "data/kd/output/{analysis}"

  # Reference strain against which group-specific variants are compared.
  # This is a biological choice: the reference defines what counts as the
  # ".0" / unsuffixed base variant.
  reference_strain: "NCIB_3610"

  # Alignment settings.
  alignment:
    substitution_matrix: "BLOSUM62"
    open_gap_score: -10
    extend_gap_score: -0.5

  # Variant classification thresholds (percentages).
  #   artefact_identity : >= this AA identity suggests the split is a
  #                       clustering artefact rather than a real variant
  #                       (Roary's own clustering threshold was 80%)
  #   min_length_ratio  : < this length ratio suggests a truncated or
  #                       fragmented gene rather than a full-length variant
  classification:
    artefact_identity: 95
    min_length_ratio: 80

# =============================================================================
# S17C - paralog / multi-copy detection
# =============================================================================
s17c:
  gbk_dir: "data/raw/gbk"
  features_file: "features_analysis_results.xlsx"   # S17A output
  output_dir: "data/kd/output/{analysis}"

  # Two copies of the same gene closer together than this are treated as a
  # candidate split gene rather than as independent paralogs.
  split_gene_distance_bp: 5000

  # Strains printed at startup as a sanity check that GBK filename mapping
  # resolved correctly. Diagnostic only; does not affect results.
  # Leave empty to skip the check.
  check_strains:
    - "NCIB_3610"
    - "PS-216"
    - "NRS6085"
    - "73"
    - "BS16045"
    - "PS-108"

# =============================================================================
# S17D - paralog annotation of curated features
# =============================================================================
s17d:
  megamatrix_67: "data/kd/inputs/67_mega_matrix_full.tsv"
  megamatrix_39: "data/kd/inputs/39_mega_matrix_full.tsv"
  corrected_file: "corrected_features_analysis.xlsx"        # S17B, reviewed
  paralog_file: "{analysis}_paralog_check_results.xlsx"     # S17C, reviewed
  output_dir: "data/kd/output/{analysis}"
  output_file: "{analysis}_corrected_features_annotated_paralogs.xlsx"

  # Strain to drop from the matrix before annotation, if any (null = none).
  excluded_strain: null

# =============================================================================
# S19 - feature enrichment and strain-centric analysis
# =============================================================================
s19:
  interaction_file: "data/strain_interactions.xlsx"
  feature_file: "data/Filtered_presence_absence_78strains.xlsx"
  output_dir: "output"

  # High-ANI subset on which the enrichment tests are run. This is the same
  # descriptive boundary used in S15A; see S15B for its validation. It is a
  # round number chosen for readability, not a fitted decision threshold.
  high_ani_threshold: 99.5

  # Inferential choices for the Fisher's exact enrichment tests.
  #   fdr_alpha        : Benjamini-Hochberg FDR cut-off for significance
  #   or_upper/or_lower: effect-size bounds for the "notable" filtered tables
  #                      (OR > upper or OR < lower)
  #   min_count        : minimum discriminating pairs carrying the feature,
  #                      applied to the filtered tables only
  #   log2or_label     : |log2 OR| above which volcano points are labelled as
  #                      directionally enriched
  enrichment:
    fdr_alpha: 0.05
    or_upper: 2.0
    or_lower: 0.5
    min_count: 5
    log2or_label: 1.0

  # Strain-centric case-control analysis. These are separate
  # from the enrichment min_count above: a strain enters the analysis only if
  # it has at least this many merging and discriminating high-ANI partners.
  strain_centric:
    min_merge: 3
    min_discriminate: 5

  plots:
    dpi: 300
    base_size: 12
    point_alpha: 0.6
    point_size: 2