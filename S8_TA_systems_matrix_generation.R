#!/usr/bin/env Rscript
# =============================================================================
# S8_TA_systems_matrix_generation.R
#
# Toxin-antitoxin (TA) systems: BLAST family classification and
# presence-absence matrix generation.
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(tidyr)
})

# --- 1. Paths -----------------------------------------------------------------
input_file <- "input/TA_pairs_results.xlsx"
output_dir <- "output"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. Load and process TA pair data -----------------------------------------
ta_pairs <- read_excel(input_file, sheet = "Split_product") %>%
  mutate(
    Ha_value       = as.numeric(Ha_value),
    E_value_BLAST  = as.numeric(E_value_BLAST),
    E_value_HMMER  = as.numeric(E_value_HMMER)
  )

cat("Strains:",          n_distinct(ta_pairs$StrainName), "\n")
cat("Total TA entries:", nrow(ta_pairs),                  "\n")

# --- 3. Create functional toxin-antitoxin pairs -------------------------------
# Each TA system consists of a toxin and antitoxin; merge into one row
# per pair with toxin always listed first.
ta_paired <- ta_pairs %>%
  mutate(Type_order = ifelse(Type == "Toxin", 1, 2)) %>%
  group_by(StrainName, Strain, TA_ID) %>%
  arrange(Type_order, .by_group = TRUE) %>%
  summarise(
    Type                    = paste(Type,                   collapse = "/"),
    Locus_tag               = paste(Locus_tag,              collapse = "/"),
    Coordinates             = paste(Coordinates,            collapse = "/"),
    Strand                  = paste(Strand,                 collapse = "/"),
    Length                  = paste(Length,                 collapse = "/"),
    Protein_ID              = paste(Protein_ID,             collapse = "/"),
    Product_gene            = paste(Product_gene,           collapse = "/"),
    Product_protein         = paste(Product_protein,        collapse = "/"),
    Product_protein_ID      = paste(Product_protein_ID,     collapse = "/"),
    TA_type                 = first(TA_type),
    `Domain(HMMER_best_hit)` = paste(`Domain(HMMER_best_hit)`, collapse = "/"),
    Blast_hit               = paste(`Blast-best_hit`,        collapse = "/"),
    Blast_family            = paste(`Blast-best_hit_family`, collapse = "/"),
    Ha_value                = paste(Ha_value,                collapse = "/"),
    E_value_BLAST           = paste(E_value_BLAST,           collapse = "/"),
    E_value_HMMER           = paste(E_value_HMMER,           collapse = "/"),
    .groups = "drop"
  )

cat("TA pairs created:", nrow(ta_paired), "\n")

# --- 4. Assign BLAST family identifiers (BF_ID) -------------------------------
# Each unique BLAST family receives a numeric identifier (BF1, BF2, ...).
# When multiple BLAST hit variants exist within the same family, they receive
# letter suffixes (BF2a, BF2b, ...).
unique_families <- ta_paired %>%
  filter(!is.na(Blast_family), Blast_family != "-/-") %>%
  distinct(Blast_family, TA_type, Product_gene, Product_protein,
           `Domain(HMMER_best_hit)`, Blast_hit) %>%
  arrange(Blast_family, Blast_hit) %>%
  group_by(Blast_family) %>%
  mutate(
    family_number = cur_group_id(),
    hit_count     = n(),
    hit_letter    = if_else(hit_count > 1, letters[row_number()], "")
  ) %>%
  ungroup() %>%
  mutate(BF_ID = paste0("BF", family_number, hit_letter)) %>%
  select(-family_number, -hit_count, -hit_letter) %>%
  select(BF_ID, everything())

cat("Unique BLAST families:",           n_distinct(unique_families$Blast_family), "\n")
cat("BLAST family variants (BF_IDs):",  nrow(unique_families),                    "\n")

write.csv(unique_families,
          file.path(output_dir, "blast_family_reference.csv"),
          row.names = FALSE)

# --- 5. Join BF_ID to TA pairs ------------------------------------------------
ta_with_bf <- ta_paired %>%
  left_join(unique_families,
            by = c("Blast_family", "TA_type", "Product_gene",
                   "Product_protein", "Domain(HMMER_best_hit)", "Blast_hit"))

cat("TA pairs with BF_ID:", sum(!is.na(ta_with_bf$BF_ID)), "of",
    nrow(ta_with_bf), "\n")

# --- 6. Generate presence-absence matrix (BLAST family variants) --------------
ta_blast_only <- ta_with_bf %>% filter(!is.na(BF_ID))

pa_matrix <- ta_blast_only %>%
  select(StrainName, BF_ID) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = BF_ID, values_from = present, values_fill = 0)

cat("Matrix dimensions:", nrow(pa_matrix), "strains x",
    ncol(pa_matrix) - 1, "BLAST family variants\n")

# --- 7. Summary statistics ----------------------------------------------------
pa_mat <- as.matrix(pa_matrix[, -1])
rownames(pa_mat) <- pa_matrix$StrainName

variant_presence <- colSums(pa_mat)
core_threshold   <- 0.95 * nrow(pa_mat)
min_threshold    <- 0.05 * nrow(pa_mat)

n_core      <- sum(variant_presence >= core_threshold)
n_accessory <- sum(variant_presence >= min_threshold &
                   variant_presence <  core_threshold)
n_rare      <- sum(variant_presence <  min_threshold)

cat("Core variants (>=95% strains):", n_core,      "\n")
cat("Accessory variants (5-95%):",    n_accessory, "\n")
cat("Rare variants (<5%):",           n_rare,      "\n")
cat("Mean TA systems per strain:",    round(mean(rowSums(pa_mat)), 2), "\n")

# --- 8. Save ------------------------------------------------------------------
output_file <- file.path(output_dir, "TA_blast_family_variants_matrix.xlsx")
write_xlsx(pa_matrix, output_file)
cat("Saved to:", output_file, "\n")
