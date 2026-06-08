#!/usr/bin/env Rscript
# =============================================================================
# S6_PHASTER_matrix_generation.R
#
# PHASTER prophage presence/absence matrix generation.
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
})

# --- 1. Paths -----------------------------------------------------------------
input_file <- "input/phaster_results.xlsx"
output_dir <- "output"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. Load PHASTER data -----------------------------------------------------
phaster_data <- read_excel(input_file, sheet = "All Phaster Results")

cat("Strains:",               length(unique(phaster_data$StrainName)), "\n")
cat("Total prophage regions:", nrow(phaster_data),                     "\n")

# --- 3. Extract main phage names ----------------------------------------------
# PHASTER reports phage hits with frequency counts
# (e.g., "PHAGE_Bacill_SPbeta_NC_001884(4)"). Keep only the phage identifier.
extract_main_phage <- function(phage_string) {
  first_phage <- strsplit(phage_string, ",")[[1]][1]
  gsub("\\([0-9]+\\)$", "", first_phage)
}

phaster_processed <- phaster_data %>%
  mutate(
    Main_Phage   = sapply(`Most Common Phage Name (all)`, extract_main_phage),
    Completeness = tolower(Completeness)
  )

cat("Unique prophage types:", length(unique(phaster_processed$Main_Phage)), "\n")

# --- 4. Completeness distribution ---------------------------------------------
completeness_summary <- phaster_processed %>%
  count(Completeness) %>%
  mutate(Percentage = round(n / sum(n) * 100, 1))

print(completeness_summary)

# --- 5. Binary presence/absence matrix (intact prophages only) ----------------
all_strains <- unique(phaster_processed$StrainName)

intact_only <- phaster_processed %>%
  filter(Completeness == "intact") %>%
  select(StrainName, Main_Phage) %>%
  distinct() %>%
  mutate(Present = 1)

matrix_intact <- intact_only %>%
  pivot_wider(
    names_from  = Main_Phage,
    values_from = Present,
    values_fill = 0
  ) %>%
  complete(StrainName = all_strains) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("StrainName")

cat("Matrix dimensions:",             nrow(matrix_intact), "strains x",
                                      ncol(matrix_intact), "prophages\n")
cat("Total intact prophage occurrences:", sum(matrix_intact),                    "\n")
cat("Mean intact prophages per strain:",  round(mean(rowSums(matrix_intact)), 2), "\n")

# --- 6. Save ------------------------------------------------------------------
output_file <- file.path(output_dir, "prophage_matrix_intact_only.xlsx")

matrix_intact %>%
  rownames_to_column("StrainName") %>%
  write_xlsx(output_file)

cat("Saved to:", output_file, "\n")
