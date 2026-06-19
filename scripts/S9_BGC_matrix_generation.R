#!/usr/bin/env Rscript
# =============================================================================
# S9_BGC_matrix_generation.R
#
# antiSMASH biosynthetic gene cluster (BGC) classification and
# presence-absence matrix generation.
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(openxlsx)
})

# --- 1. Paths -----------------------------------------------------------------
antismash_file <- "input/antismash_results.xlsx"
output_dir     <- "output"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. Load input data -------------------------------------------------------
antismash_data <- read_excel(antismash_file, sheet = "Sheet 1")
chrom_lengths  <- read_excel(antismash_file, sheet = "chromosome_length")

cat("Strains:",           n_distinct(antismash_data$StrainName), "\n")
cat("Total BGC regions:", nrow(antismash_data),                  "\n")

# --- 3. Clean coordinates and validate ----------------------------------------
antismash_data <- antismash_data %>%
  mutate(
    From          = as.numeric(str_remove_all(as.character(From),          "[^0-9]")),
    To            = as.numeric(str_remove_all(as.character(To),            "[^0-9]")),
    Region_length = as.numeric(str_remove_all(as.character(Region_length), "[^0-9]"))
  ) %>%
  filter(!is.na(From), !is.na(To), !is.na(Region_length),
         From < To, Region_length > 0)

cat("Valid BGC regions after cleaning:", nrow(antismash_data), "\n")

# --- 4. Classify BGCs into three tiers ----------------------------------------
antismash_processed <- antismash_data %>%
  mutate(
    Similarity_percent = as.numeric(gsub("%", "", Similarity)),
    BGC_Classification = case_when(
      !is.na(Most_similar_known_cluster) & Similarity_percent >= 70 ~ "Known",
      !is.na(Most_similar_known_cluster) & Similarity_percent >= 30 ~ "LowConf",
      TRUE                                                          ~ "Novel"
    )
  )

classification_summary <- antismash_processed %>%
  count(BGC_Classification) %>%
  mutate(Percentage = round(n / sum(n) * 100, 1))

print(classification_summary)

# --- 5. Name known and low-confidence BGCs ------------------------------------
known_lowconf <- antismash_processed %>%
  filter(BGC_Classification %in% c("Known", "LowConf")) %>%
  mutate(
    BGC_Name = case_when(
      BGC_Classification == "Known" ~
        gsub("[^A-Za-z0-9_]", "_", Most_similar_known_cluster),
      BGC_Classification == "LowConf" ~
        paste0("LowConf_", gsub("[^A-Za-z0-9_]", "_", Most_similar_known_cluster))
    )
  )

# --- 6. Cluster novel BGCs by type, size, and chromosomal position ------------
# Novel BGCs are grouped using three orthogonal features:
#   - Biosynthetic pathway type (from antiSMASH `Type` field)
#   - Region size in 10-kb bins
#   - Relative chromosomal position in 20% bins
# Groups with >= 3 occurrences retain the full descriptor; smaller groups
# receive a simplified identifier.
novel_clustered <- antismash_processed %>%
  filter(BGC_Classification == "Novel") %>%
  left_join(
    chrom_lengths %>% select(Header_ID, Chromosome_Length),
    by = c("Strain" = "Header_ID")
  ) %>%
  filter(!is.na(Chromosome_Length)) %>%
  mutate(
    Region_midpoint           = (From + To) / 2,
    Relative_position_percent = (Region_midpoint / Chromosome_Length) * 100,
    Size_bin                  = pmax(round(Region_length / 10000) * 10000, 10000),
    Position_bin              = pmax(0, pmin(100, round(Relative_position_percent / 20) * 20)),
    Clean_type                = gsub("[^A-Za-z0-9]", "_", Type),
    Cluster_group             = paste0(Clean_type, "_", Size_bin / 1000, "kb_",
                                       Position_bin, "pct")
  ) %>%
  group_by(Cluster_group) %>%
  mutate(Group_size = n(), Group_ID = cur_group_id()) %>%
  ungroup() %>%
  mutate(
    BGC_Name = case_when(
      Group_size >= 3 ~ paste0("Novel_", Cluster_group),
      TRUE            ~ paste0("Novel_", Clean_type, "_", Group_ID)
    )
  )

# Novel BGCs that could not be positioned (missing chromosome length)
novel_unpositioned <- antismash_processed %>%
  filter(BGC_Classification == "Novel") %>%
  anti_join(novel_clustered, by = c("StrainName", "Region")) %>%
  mutate(
    BGC_Name = paste0("Novel_", gsub("[^A-Za-z0-9]", "_", Type),
                      "_unpositioned_", row_number())
  )

cat("Novel BGC groups after clustering:", n_distinct(novel_clustered$BGC_Name), "\n")
cat("Unpositioned novel BGCs:",           nrow(novel_unpositioned),              "\n")

# --- 7. Combine all BGCs ------------------------------------------------------
antismash_named <- bind_rows(
  known_lowconf      %>% select(all_of(names(antismash_processed)), BGC_Name),
  novel_clustered    %>% select(all_of(names(antismash_processed)), BGC_Name),
  novel_unpositioned %>% select(all_of(names(antismash_processed)), BGC_Name)
)

final_summary <- antismash_named %>%
  group_by(BGC_Classification) %>%
  summarise(Unique_BGCs = n_distinct(BGC_Name),
            Occurrences = n(),
            .groups     = "drop")

print(final_summary)

# --- 8. Cluster validation ----------------------------------------------------
# Assess consistency of BGCs sharing the same name across strains using
# coefficients of variation (CV) for region length and chromosomal position.
# Clusters with CV > 30% are flagged for manual review.
antismash_positioned <- antismash_named %>%
  left_join(
    chrom_lengths %>% select(Header_ID, Chromosome_Length),
    by = c("Strain" = "Header_ID")
  ) %>%
  filter(!is.na(Chromosome_Length)) %>%
  mutate(
    Region_midpoint           = (From + To) / 2,
    Relative_position_percent = (Region_midpoint / Chromosome_Length) * 100
  )

cluster_validation <- antismash_positioned %>%
  group_by(BGC_Name, BGC_Classification) %>%
  summarise(
    N_strains       = n_distinct(StrainName),
    Count           = n(),
    Mean_length     = mean(Region_length,              na.rm = TRUE),
    SD_length       = sd(Region_length,                na.rm = TRUE),
    Min_length      = min(Region_length,               na.rm = TRUE),
    Max_length      = max(Region_length,               na.rm = TRUE),
    CV_length       = (SD_length / Mean_length) * 100,
    Mean_position   = mean(Relative_position_percent,  na.rm = TRUE),
    SD_position     = sd(Relative_position_percent,    na.rm = TRUE),
    CV_position     = (SD_position / Mean_position) * 100,
    Mean_similarity = mean(Similarity_percent,         na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  mutate(
    Length_flag = case_when(
      CV_length > 30 & N_strains > 1 ~ "High variation",
      is.na(CV_length) & N_strains > 1 ~ "Check manually",
      N_strains == 1                 ~ "Single strain",
      TRUE                           ~ "OK"
    )
  ) %>%
  arrange(desc(N_strains))

flagged <- cluster_validation %>%
  filter(Length_flag %in% c("High variation", "Check manually"), N_strains > 1)

cat("Total unique BGCs:",       nrow(cluster_validation),              "\n")
cat("BGCs in multiple strains:", sum(cluster_validation$N_strains > 1), "\n")
cat("Flagged (CV > 30%):",       nrow(flagged),                         "\n")

write_xlsx(cluster_validation, file.path(output_dir, "BGC_cluster_validation.xlsx"))

# --- 9. Generate presence-absence matrix --------------------------------------
bgc_matrix <- antismash_named %>%
  select(StrainName, BGC_Name) %>%
  distinct() %>%
  mutate(Present = 1) %>%
  pivot_wider(names_from = BGC_Name, values_from = Present, values_fill = 0)

bgc_classifications <- antismash_named %>%
  select(BGC_Name, BGC_Classification) %>%
  distinct()

cat("Matrix:", nrow(bgc_matrix), "strains x", ncol(bgc_matrix) - 1, "BGCs\n")
cat("  Known:",   sum(bgc_classifications$BGC_Classification == "Known"),   "\n")
cat("  LowConf:", sum(bgc_classifications$BGC_Classification == "LowConf"), "\n")
cat("  Novel:",   sum(bgc_classifications$BGC_Classification == "Novel"),   "\n")

bgc_matrix_export <- bgc_matrix %>% column_to_rownames("StrainName")

# --- 10. Save outputs ---------------------------------------------------------
# Colour-coded matrix: header cells are tinted by BGC classification.
wb <- createWorkbook()
addWorksheet(wb, "BGC_Matrix")
writeData(wb, "BGC_Matrix", bgc_matrix_export, rowNames = TRUE)

col_positions <- data.frame(
  BGC_Name = colnames(bgc_matrix_export),
  Col_num  = 2:(ncol(bgc_matrix_export) + 1)
) %>%
  left_join(bgc_classifications, by = "BGC_Name")

for (cls in c("Known", "LowConf", "Novel")) {
  cols <- col_positions %>% filter(BGC_Classification == cls) %>% pull(Col_num)
  fill_color <- switch(cls,
                       Known   = "#D4EDDA",
                       LowConf = "#FFF3CD",
                       Novel   = "#FFE5CC")
  if (length(cols) > 0) {
    addStyle(wb, "BGC_Matrix",
             createStyle(fgFill = fill_color),
             rows = 1, cols = cols,
             gridExpand = TRUE, stack = TRUE)
  }
}

addStyle(wb, "BGC_Matrix", createStyle(textDecoration = "bold"),
         rows = 1, cols = 1:(ncol(bgc_matrix_export) + 1),
         gridExpand = TRUE, stack = TRUE)

saveWorkbook(wb, file.path(output_dir, "BGC_matrix.xlsx"), overwrite = TRUE)

# Detailed annotations table
bgc_details <- antismash_named %>%
  select(StrainName, Strain, Region, BGC_Name, BGC_Classification, Type,
         Most_similar_known_cluster, Similarity_percent,
         Region_length, From, To) %>%
  arrange(StrainName, Region)

write_xlsx(bgc_details, file.path(output_dir, "BGC_detailed_annotations.xlsx"))

cat("Saved: BGC_matrix.xlsx, BGC_cluster_validation.xlsx, BGC_detailed_annotations.xlsx\n")
