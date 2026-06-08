#!/usr/bin/env Rscript
# =============================================================================
# S17_feature_enrichment_and_strain_centric_analysis.R
#
# Feature enrichment and strain-centric case-control analysis of kin
# discrimination determinants (Methods 5.11 and 5.12).
#
# Analyses performed:
#   (1) Feature enrichment across 399 high-ANI (fastANI >= 99.5%) pairwise
#       interactions using Fisher's exact tests for
#         (a) shared presence (both strains carry the feature) and
#         (b) feature mismatch (XOR: exactly one strain carries it),
#       with Benjamini-Hochberg FDR correction.
#       Produces Tables S26, S27A-B, S28, S29A-B and the volcano plots
#       shown in Supplementary Figure S6.
#   (2) Strain-centric case-control analysis for focal strains with
#       >=3 merging and >=5 discriminating high-ANI partners, classifying
#       features into four presence/absence categories.
#       Produces Tables S30 and S31.
#
# Inputs:
#   data/strain_interactions.xlsx            -- pairwise swarming + fastANI
#                                               data (850 pairs).
#   data/Filtered_presence_absence_78strains.xlsx
#                                            -- multi-sheet feature matrix
#                                               (genes, prophages, TA systems,
#                                               TFs, metabolic reactions, BGCs).
#
# Threshold summaries, Pearson correlation, and nonconforming-pair tabulation
# are produced separately by S13_strain_compatibility_analysis.R.
#
# Usage : Rscript S17_feature_enrichment_and_strain_centric_analysis.R
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(kableExtra)
  library(ggrepel)
})

# --- Output directory ---------------------------------------------------------
run_output_dir <- "output"
dir.create(run_output_dir,                       recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 12))

# --- Paths (edit as needed) ---------------------------------------------------
interaction_file <- "data/strain_interactions.xlsx"
feature_file     <- "data/Filtered_presence_absence_78strains.xlsx"

# Helper: rename strain identifiers via a lookup vector
rename_strains <- function(x, key) ifelse(x %in% names(key), key[x], x)

# --- 1. Data loading and preparation ------------------------------------------

# --- Load interactions --------------------------------------------------------
data <- read_excel(interaction_file)

required_cols <- c("S1", "S2", "Phenotype_interaction",
                   "Group_S1", "Group_S2", "fastANI")
missing_cols  <- setdiff(required_cols, names(data))
if (length(missing_cols)) stop("Missing columns: ",
                               paste(missing_cols, collapse = ", "))

cat("Pairs loaded:", nrow(data), "\n")

# --- Transform phenotypes -----------------------------------------------------
convert_phenotype <- function(p) {
  case_when(
    str_detect(p, "Merging") | p == "1"    ~ 1.0,
    str_detect(p, "Contact") | p == "0.5"  ~ 0.5,
    str_detect(p, "Thick")   | p == "0.25" ~ 0.25,
    str_detect(p, "Lysis")   | p == "0"    ~ 0.0,
    TRUE ~ as.numeric(str_extract(p, "[0-9.]+"))
  )
}

data_clean <- data |>
  mutate(
    S1                 = as.character(S1),
    S2                 = as.character(S2),
    phenotype_numeric  = convert_phenotype(Phenotype_interaction),
    same_group         = Group_S1 == Group_S2,
    fastANI_threshold  = fastANI >= 99.5,
    observed_kin       = phenotype_numeric >= 0.5,
    phenotype_category = case_when(
      phenotype_numeric == 1.0  ~ "Merging",
      phenotype_numeric == 0.5  ~ "Contact",
      phenotype_numeric == 0.25 ~ "Thick boundary",
      phenotype_numeric == 0.0  ~ "Lysis",
      TRUE                      ~ "Unknown"
    )
  )

# --- Harmonize strain identifiers ---------------------------------------------
# Feature-style (left) -> pair-style (right). Feature matrix columns use
# compact names (PS216); interaction data uses punctuated names (PS-216).
strain_rename_key <- c(
  "NCIB_3610" = "3610", "BS16045" = "Bs16045",
  "MB8_B1" = "MB8B1",  "MB8_B10" = "MB8B10", "MB8_B7" = "MB8B7",
  "MB9_B1" = "MB9B1",  "MB9_B4"  = "MB9B4",  "MB9_B6" = "MB9B6",
  "P8_B1"  = "P8B1",   "P8_B3"   = "P8B3",   "P9_B1"  = "P9B1",
  "PS-11" = "PS11",   "PS-13" = "PS13",   "PS-14" = "PS14",   "PS-15" = "PS15",
  "PS-18" = "PS18",   "PS-20" = "PS20",   "PS-24" = "PS24",   "PS-25" = "PS25",
  "PS-30" = "PS30",   "PS-31" = "PS31",   "PS-51" = "PS51",   "PS-52" = "PS52",
  "PS-53" = "PS53",   "PS-54" = "PS54",   "PS-55" = "PS55",   "PS-64" = "PS64",
  "PS-65" = "PS65",   "PS-68" = "PS68",   "PS-93" = "PS93",   "PS-95" = "PS95",
  "PS-96" = "PS96",   "PS-108" = "PS108", "PS-109" = "PS109", "PS-119" = "PS119",
  "PS-130" = "PS130", "PS-131" = "PS131", "PS-149" = "PS149", "PS-160" = "PS160",
  "PS-168" = "PS168", "PS-194" = "PS194", "PS-196" = "PS196", "PS-209" = "PS209",
  "PS-210" = "PS210", "PS-216" = "PS216", "PS-217" = "PS217", "PS-218" = "PS218",
  "PS-233" = "PS233", "PS-237" = "PS237", "PS-261" = "PS261", "PS-263" = "PS263",
  "KF24" = "KF24", "73" = "73",
  "NRS6085" = "NRS6085", "NRS6099" = "NRS6099", "NRS6103" = "NRS6103",
  "NRS6105" = "NRS6105", "NRS6107" = "NRS6107", "NRS6116" = "NRS6116",
  "NRS6118" = "NRS6118", "NRS6121" = "NRS6121", "NRS6127" = "NRS6127",
  "NRS6128" = "NRS6128", "NRS6132" = "NRS6132", "NRS6145" = "NRS6145",
  "NRS6153" = "NRS6153", "NRS6160" = "NRS6160", "NRS6181" = "NRS6181",
  "NRS6183" = "NRS6183", "NRS6186" = "NRS6186", "NRS6187" = "NRS6187",
  "NRS6190" = "NRS6190", "NRS6202" = "NRS6202",
  "RO-F-3" = "ROF3", "RS-D-2" = "RSD2", "RO-DD-2" = "RODD2",
  "RO-A-4" = "ROA4", "RO-FF-1" = "ROFF1"
)

data_clean <- data_clean |>
  mutate(S1 = rename_strains(S1, strain_rename_key),
         S2 = rename_strains(S2, strain_rename_key))

# --- Load feature matrix ------------------------------------------------------
sheet_names <- excel_sheets(feature_file)

mega_matrix <- lapply(sheet_names, read_excel, path = feature_file)
names(mega_matrix) <- sheet_names

# Drop non-matrix sheets and rename strain columns to pair-style
mega_matrix_renamed <- mega_matrix[!names(mega_matrix) %in% c("Reduction_Summary")]
mega_matrix_renamed <- lapply(mega_matrix_renamed, function(df) {
  colnames(df)[-1] <- rename_strains(colnames(df)[-1], strain_rename_key)
  df
})

features_loaded <- length(mega_matrix_renamed) > 0
cat("Feature classes loaded:", length(mega_matrix_renamed),              "\n")
cat("Total features:",         sum(sapply(mega_matrix_renamed, nrow)),   "\n")

# --- 2. Feature enrichment (high-ANI pairs) -----------------------------------

# --- 2.1 Define target and background set -------------------------------------
data_annotated <- data_clean |>
  mutate(HighANI_Discriminating = fastANI >= 99.5 & !observed_kin)

data_highANI <- data_annotated |> filter(fastANI >= 99.5)

cat("High-ANI pairs (fastANI >= 99.5%):       ",
    nrow(data_highANI),                       "\n")
cat("High-ANI discriminating (nonconforming): ",
    sum(data_highANI$HighANI_Discriminating), "\n")

# --- 2.2 Test functions -------------------------------------------------------

# Fisher's exact test on feature *shared presence* (both strains carry the
# feature) vs HighANI_Discriminating status.
test_feature_enrichment <- function(feature_vec, interaction_df) {
  both_feature <- interaction_df$S1 %in% feature_vec &
                  interaction_df$S2 %in% feature_vec
  tbl <- table(HighANI_Discriminating = interaction_df$HighANI_Discriminating,
               Both_Feature           = both_feature)
  if (min(dim(tbl)) < 2) return(NULL)
  ft <- fisher.test(tbl)
  tibble(
    OR                = unname(ft$estimate),
    CI_low            = ft$conf.int[1],
    CI_high           = ft$conf.int[2],
    p_value           = ft$p.value,
    n_highANI_disc    = sum(interaction_df$HighANI_Discriminating & both_feature),
    pct_highANI_disc  = mean(both_feature[ interaction_df$HighANI_Discriminating]) * 100,
    pct_background    = mean(both_feature[!interaction_df$HighANI_Discriminating]) * 100,
    n_feature_strains = length(feature_vec)
  )
}

# Fisher's exact test on feature *mismatch* (XOR: exactly one strain
# carries the feature).
compute_feature_mismatch <- function(feature_vec, interaction_df) {
  s1_has   <- interaction_df$S1 %in% feature_vec
  s2_has   <- interaction_df$S2 %in% feature_vec
  mismatch <- xor(s1_has, s2_has)
  tbl <- table(HighANI_Discriminating = interaction_df$HighANI_Discriminating,
               Mismatch               = mismatch)
  if (min(dim(tbl)) < 2) return(NULL)
  ft <- fisher.test(tbl)
  tibble(
    OR                = unname(ft$estimate),
    CI_low            = ft$conf.int[1],
    CI_high           = ft$conf.int[2],
    p_value           = ft$p.value,
    n_mismatch_disc   = sum(interaction_df$HighANI_Discriminating & mismatch),
    pct_mismatch_disc = mean(mismatch[ interaction_df$HighANI_Discriminating]) * 100,
    pct_mismatch_bg   = mean(mismatch[!interaction_df$HighANI_Discriminating]) * 100,
    n_feature_strains = length(feature_vec)
  )
}

# --- 2.3 Shared-presence enrichment -------------------------------------------
if (features_loaded) {

  feature_results <- lapply(names(mega_matrix_renamed), function(sheet) {
    df          <- mega_matrix_renamed[[sheet]]
    feature_col <- colnames(df)[1]
    strain_cols <- colnames(df)[-1]

    long_df <- df |>
      pivot_longer(all_of(strain_cols),
                   names_to  = "Strain",
                   values_to = "Present") |>
      filter(Present > 0)

    features <- unique(long_df[[feature_col]])

    bind_rows(lapply(features, function(feat) {
      strains_with_feature <- long_df |>
        filter(.data[[feature_col]] == feat) |>
        pull(Strain) |>
        unique()
      out <- test_feature_enrichment(strains_with_feature, data_highANI)
      if (is.null(out)) return(NULL)
      out |> mutate(Feature = feat, Class = sheet)
    }))
  })

  feature_results_all <- bind_rows(feature_results) |>
    mutate(FDR = p.adjust(p_value, method = "fdr")) |>
    arrange(FDR, desc(OR))

  cat("Shared-presence features tested:",
      nrow(feature_results_all),                  "\n")
  cat("Significant (FDR < 0.05):       ",
      sum(feature_results_all$FDR < 0.05),        "\n")

  # Stringent candidates: |log2OR| > 1 AND n >= 5
  stringent_shared <- feature_results_all |>
    filter(FDR < 0.05, (OR > 2 | OR < 0.5), n_highANI_disc >= 5,
           is.finite(OR), OR > 0) |>
    mutate(
      logOR     = log2(OR),
      abs_logOR = abs(logOR),
      Effect    = ifelse(OR > 1, "up_Disc", "up_Kin"),
      OR_CI     = sprintf("%.2f [%.2f-%.2f]", OR, CI_low, CI_high)
    )

  cat("Stringent shared candidates (|log2OR|>1, n>=5):",
      nrow(stringent_shared), "\n")

  # ---- Exports ---------------------------------------------------------------
  # Table S26: all significant shared-presence features (FDR < 0.05)
  write_csv(feature_results_all |> filter(FDR < 0.05),
            file.path(run_output_dir, "TableS26_shared_significant.csv"))

  # Table S27A: top 20 stringent candidates (FDR-first)
  write_csv(stringent_shared |>
              arrange(FDR, desc(abs_logOR), desc(n_highANI_disc)) |>
              slice_head(n = 20),
            file.path(run_output_dir, "TableS27A_shared_top20.csv"))

  # Table S27B: top 5 stringent candidates per feature class
  write_csv(stringent_shared |>
              group_by(Class) |>
              arrange(FDR, desc(abs_logOR), desc(n_highANI_disc),
                      .by_group = TRUE) |>
              slice_head(n = 5) |>
              ungroup(),
            file.path(run_output_dir, "TableS27B_shared_top5_per_class.csv"))

  # ---- Volcano plot ----------------------------------------------------------
  volcano_shared <- feature_results_all |>
    mutate(
      log2OR        = log2(OR),
      neg_log10_FDR = -log10(FDR),
      Significance  = case_when(
        FDR < 0.05 & log2OR >  1 ~ "Enriched in discrimination",
        FDR < 0.05 & log2OR < -1 ~ "Enriched in kin behaviour",
        FDR < 0.05                ~ "Significant (small effect)",
        TRUE                      ~ "Not significant"
      ),
      label = ifelse(FDR < 0.01 & abs(log2OR) > 2.5,
                     ifelse(grepl("^hypothetical", Feature, ignore.case = TRUE),
                            paste0("hypo_", sub(".*_", "", Feature)),
                            Feature),
                     NA)
    )

  p_volcano_shared <- ggplot(volcano_shared, aes(log2OR, neg_log10_FDR)) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", color = "gray40") +
    geom_text_repel(aes(label = label), size = 2.5,
                    max.overlaps = 15, box.padding = 0.3) +
    scale_color_manual(values = c(
      "Enriched in discrimination" = "#E41A1C",
      "Enriched in kin behaviour"  = "#377EB8",
      "Significant (small effect)" = "#984EA3",
      "Not significant"            = "gray70"
    )) +
    labs(
      x     = expression(log[2](Odds~Ratio)),
      y     = expression(-log[10](FDR)),
      color = "Category"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  print(p_volcano_shared)
  ggsave(file.path(run_output_dir, "volcano_shared.svg"),
         p_volcano_shared, width = 8, height = 6, units = "in")

  # =============================================================================
  # 2.4 Mismatch enrichment
  # =============================================================================
  mismatch_results_by_class <- list()

  for (sheet in names(mega_matrix_renamed)) {
    df          <- mega_matrix_renamed[[sheet]]
    feature_col <- colnames(df)[1]
    strain_cols <- colnames(df)[-c(1, which(colnames(df) %in% c("Strain_raw", "Strain")))]

    class_results <- list()
    for (i in seq_len(nrow(df))) {
      feature_name         <- df[[feature_col]][i]
      strains_with_feature <- strain_cols[as.numeric(df[i, strain_cols]) == 1]
      res <- compute_feature_mismatch(strains_with_feature, data_highANI)
      if (!is.null(res)) class_results[[feature_name]] <- res
    }

    mismatch_results_by_class[[sheet]] <-
      bind_rows(class_results, .id = "Feature") |>
      mutate(Feature_Class = sheet)
  }

  mismatch_results_all <- bind_rows(mismatch_results_by_class) |>
    mutate(FDR = p.adjust(p_value, method = "BH"))

  cat("Mismatch features tested:   ",
      nrow(mismatch_results_all),                 "\n")
  cat("Significant (FDR < 0.05):   ",
      sum(mismatch_results_all$FDR < 0.05),       "\n")

  # Stringent mismatch candidates
  stringent_mismatch <- mismatch_results_all |>
    filter(FDR < 0.05, (OR > 2 | OR < 0.5), n_mismatch_disc >= 5,
           is.finite(OR), OR > 0) |>
    mutate(
      logOR     = log2(OR),
      abs_logOR = abs(logOR),
      Effect    = ifelse(OR > 1, "up_Disc", "up_Kin"),
      OR_CI     = sprintf("%.2f [%.2f-%.2f]", OR, CI_low, CI_high)
    )

  cat("Stringent mismatch candidates (|log2OR|>1, n>=5):",
      nrow(stringent_mismatch), "\n")

  # ---- Exports ---------------------------------------------------------------
  # Table S28: all significant mismatch features (FDR < 0.05)
  write_csv(mismatch_results_all |> filter(FDR < 0.05),
            file.path(run_output_dir, "TableS28_mismatch_significant.csv"))

  # Table S29A: top 20 stringent mismatch candidates (FDR-first)
  write_csv(stringent_mismatch |>
              arrange(FDR, desc(abs_logOR), desc(n_mismatch_disc)) |>
              slice_head(n = 20),
            file.path(run_output_dir, "TableS29A_mismatch_top20.csv"))

  # Table S29B: top 5 stringent mismatch candidates per feature class
  write_csv(stringent_mismatch |>
              group_by(Feature_Class) |>
              arrange(FDR, desc(abs_logOR), desc(n_mismatch_disc),
                      .by_group = TRUE) |>
              slice_head(n = 5) |>
              ungroup(),
            file.path(run_output_dir, "TableS29B_mismatch_top5_per_class.csv"))

  # ---- Volcano plot ----------------------------------------------------------
  x_limit <- 4.5
  y_limit <- 14

  volcano_mismatch <- mismatch_results_all |>
    mutate(
      log2OR             = log2(OR),
      neg_log10_FDR      = -log10(FDR),
      log2OR_plot        = pmin(pmax(log2OR, -x_limit), x_limit),
      neg_log10_FDR_plot = pmin(neg_log10_FDR, y_limit),
      clipped            = log2OR >  x_limit | log2OR < -x_limit |
                           neg_log10_FDR > y_limit,
      Significance       = case_when(
        FDR < 0.05 & log2OR >  1 ~ "Mismatch promotes discrimination",
        FDR < 0.05 & log2OR < -1 ~ "Mismatch promotes kin behaviour",
        FDR < 0.05                ~ "Significant (small effect)",
        TRUE                      ~ "Not significant"
      ),
      label = ifelse(FDR < 0.01 & abs(log2OR) > 3,
                     ifelse(grepl("^hypothetical", Feature, ignore.case = TRUE),
                            paste0("hypo_", sub(".*_", "", Feature)),
                            Feature),
                     NA)
    )

  p_volcano_mismatch <- ggplot(volcano_mismatch,
                               aes(log2OR_plot, neg_log10_FDR_plot)) +
    geom_point(aes(color = Significance, shape = clipped),
               alpha = 0.6, size = 2) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                       guide = "none") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", color = "gray40") +
    geom_text_repel(aes(label = label), size = 2.5,
                    max.overlaps = 20, box.padding = 0.3) +
    scale_color_manual(values = c(
      "Mismatch promotes discrimination" = "#E41A1C",
      "Mismatch promotes kin behaviour"  = "#377EB8",
      "Significant (small effect)"       = "#984EA3",
      "Not significant"                  = "gray70"
    )) +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = seq(-4, 4, by = 1),
      labels = function(x) ifelse(x ==  x_limit, paste0(">=", x),
                           ifelse(x == -x_limit, paste0("<=", x), x))
    ) +
    scale_y_continuous(limits = c(0, y_limit)) +
    labs(
      x     = expression(log[2](Odds~Ratio)),
      y     = expression(-log[10](FDR)),
      color = "Category"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  print(p_volcano_mismatch)
  ggsave(file.path(run_output_dir, "volcano_mismatch.svg"),
         p_volcano_mismatch, width = 9, height = 7, units = "in")
}

# --- 3. Strain-centric case-control analysis ----------------------------------

# --- 3.1 Identify focal strains -----------------------------------------------
# Focal strains: >=3 merging AND >=5 discriminating high-ANI partners.
high_ANI_all <- data_clean |> filter(fastANI >= 99.5)

strain_partners <- bind_rows(
  high_ANI_all |> select(Strain = S1, Partner = S2,
                         fastANI, observed_kin, phenotype_category),
  high_ANI_all |> select(Strain = S2, Partner = S1,
                         fastANI, observed_kin, phenotype_category)
) |>
  group_by(Strain) |>
  summarise(
    n_partners     = n(),
    n_merge        = sum( observed_kin),
    n_discriminate = sum(!observed_kin),
    merge_partners = paste(Partner[ observed_kin], collapse = ", "),
    disc_partners  = paste(Partner[!observed_kin], collapse = ", "),
    .groups = "drop"
  )

focal_strains_df <- strain_partners |>
  filter(n_merge >= 3, n_discriminate >= 5) |>
  arrange(desc(n_discriminate))

cat("Focal strains selected:", nrow(focal_strains_df), "\n")
print(focal_strains_df |>
        select(Strain, n_partners, n_merge, n_discriminate))

# --- 3.2 Per-focal-strain feature classification ------------------------------
# Category 1 (strict):   focal HAS;  <=20% disc, >=60% merge
# Category 2 (strict):   focal LACKS; >=80% disc, <=40% merge
# Category 3 (relaxed):  focal HAS;  <=30% disc, fewer disc than merge
# Category 4 (relaxed):  focal LACKS; >=70% disc, fewer merge than disc

# All features that appear anywhere in the matrix
all_genes_in_dataset <- unique(unlist(
  lapply(mega_matrix_renamed, function(df) df[[colnames(df)[1]]])
))

# Presence lookup: which features does a given strain carry?
get_strain_genes <- function(strain, matrices) {
  unique(unlist(lapply(matrices, function(df) {
    feat_col <- colnames(df)[1]
    if (!strain %in% colnames(df)) return(character())
    df[[feat_col]][as.numeric(df[[strain]]) == 1]
  })))
}

results_list <- list()

for (i in seq_len(nrow(focal_strains_df))) {
  focal          <- focal_strains_df$Strain[i]
  disc_partners  <- str_split(focal_strains_df$disc_partners[i],  ", ")[[1]]
  merge_partners <- str_split(focal_strains_df$merge_partners[i], ", ")[[1]]

  focal_genes <- get_strain_genes(focal, mega_matrix_renamed)
  disc_partner_genes <- setNames(
    lapply(disc_partners,  get_strain_genes,
           matrices = mega_matrix_renamed),
    disc_partners
  )
  merge_partner_genes <- setNames(
    lapply(merge_partners, get_strain_genes,
           matrices = mega_matrix_renamed),
    merge_partners
  )

  gene_stats <- data.frame(Gene = all_genes_in_dataset) |>
    rowwise() |>
    mutate(
      focal_has     = Gene %in% focal_genes,
      n_disc_has    = sum(sapply(disc_partner_genes,
                                 function(g) Gene %in% g)),
      n_merge_has   = sum(sapply(merge_partner_genes,
                                 function(g) Gene %in% g)),
      pct_disc_has  = round(n_disc_has  / length(disc_partners)  * 100, 1),
      pct_merge_has = round(n_merge_has / length(merge_partners) * 100, 1)
    ) |>
    ungroup()

  results_list[[focal]] <- list(
    focal   = focal,
    n_disc  = length(disc_partners),
    n_merge = length(merge_partners),
    focal_has_disc_lack   = gene_stats |>
      filter(focal_has,  pct_disc_has <= 20, pct_merge_has >= 60) |>
      arrange(pct_disc_has, desc(pct_merge_has)),
    focal_lacks_disc_have = gene_stats |>
      filter(!focal_has, pct_disc_has >= 80, pct_merge_has <= 40) |>
      arrange(desc(pct_disc_has), pct_merge_has),
    focal_has_mismatch    = gene_stats |>
      filter(focal_has,  pct_disc_has <= 30, n_disc_has  < n_merge_has) |>
      arrange(pct_disc_has),
    focal_lacks_mismatch  = gene_stats |>
      filter(!focal_has, pct_disc_has >= 70, n_merge_has < n_disc_has) |>
      arrange(desc(pct_disc_has))
  )
}

# --- 3.3 Compile candidates and identify recurrent features -------------------
all_candidates <- bind_rows(lapply(names(results_list), function(focal) {
  res <- results_list[[focal]]
  bind_rows(
    if (nrow(res$focal_has_disc_lack))   res$focal_has_disc_lack   |>
        mutate(Focal = focal,
               Category = "Cat1: Focal HAS (strict)",
               Direction = "Focal HAS"),
    if (nrow(res$focal_lacks_disc_have)) res$focal_lacks_disc_have |>
        mutate(Focal = focal,
               Category = "Cat2: Focal LACKS (strict)",
               Direction = "Focal LACKS"),
    if (nrow(res$focal_has_mismatch))    res$focal_has_mismatch    |>
        mutate(Focal = focal,
               Category = "Cat3: Focal HAS (relaxed)",
               Direction = "Focal HAS"),
    if (nrow(res$focal_lacks_mismatch))  res$focal_lacks_mismatch  |>
        mutate(Focal = focal,
               Category = "Cat4: Focal LACKS (relaxed)",
               Direction = "Focal LACKS")
  )
}))

cat("Total feature-focal combinations:",
    nrow(all_candidates),                "\n")
cat("Unique candidate features:       ",
    length(unique(all_candidates$Gene)), "\n")

gene_freq <- all_candidates |>
  distinct(Gene, Focal, Direction) |>
  count(Gene, Direction, name = "N_Focal_Strains") |>
  arrange(desc(N_Focal_Strains))

recurrent <- gene_freq |> filter(N_Focal_Strains >= 2)

cat("Recurrent (>=2 focal strains):\n")
cat("  Focal HAS:  ", sum(recurrent$Direction == "Focal HAS"),   "\n")
cat("  Focal LACKS:", sum(recurrent$Direction == "Focal LACKS"), "\n")

# --- Exports ------------------------------------------------------------------
# Table S30: per-focal-strain candidate counts and partner composition
focal_summary <- focal_strains_df |>
  select(Strain, n_partners, n_merge, n_discriminate) |>
  left_join(
    all_candidates |>
      count(Focal, Category) |>
      pivot_wider(names_from = Category, values_from = n, values_fill = 0),
    by = c("Strain" = "Focal")
  )
write_csv(focal_summary,
          file.path(run_output_dir, "TableS30_focal_strain_summary.csv"))

# Full candidate table (all focal-feature combinations)
write_csv(all_candidates,
          file.path(run_output_dir, "strain_centric_all_candidates.csv"))

# Table S31: recurrent candidates (>=2 focal strain contexts) with
# supporting focal strain lists
recurrent_with_strains <- recurrent |>
  left_join(
    all_candidates |>
      group_by(Gene, Direction) |>
      summarise(Focal_Strains = paste(unique(Focal), collapse = ", "),
                .groups = "drop"),
    by = c("Gene", "Direction")
  )
write_csv(recurrent_with_strains,
          file.path(run_output_dir, "TableS31_recurrent_candidates.csv"))
