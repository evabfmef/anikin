#!/usr/bin/env Rscript
# =============================================================================
# S11B_TF_matrix_generation.R
#
# Transcription factor (TF) presence-absence matrix generation.
#
# Processes PredicTF output files to construct a binary presence-absence
# matrix of transcription factor (TF) families across strains.
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(writexl)
})

# --- 1. Paths -----------------------------------------------------------------
files_directory <- "input/collected_mapping.potential.TF"
output_dir      <- "output"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. Load PredicTF output files --------------------------------------------
file_list <- list.files(files_directory, full.names = TRUE, recursive = FALSE)
file_list <- file_list[!file.info(file_list)$isdir]

cat("Total files found:", length(file_list), "\n")

# --- 3. Extract unique TF families per genome ---------------------------------
accessions_dict <- list()

for (file_path in file_list) {
  lines      <- readLines(file_path)
  data_lines <- lines[!grepl("^#", lines) & nchar(lines) > 0]

  tf_categories <- sapply(strsplit(data_lines, "\\s+"), `[`, 1)
  unique_tfs    <- unique(tf_categories)

  genome_name <- strsplit(basename(file_path), "\\.")[[1]][1]
  accessions_dict[[genome_name]] <- unique_tfs
}

all_TFs <- unique(unlist(accessions_dict))

cat("Genomes processed:",   length(accessions_dict), "\n")
cat("Unique TF families:",  length(all_TFs),         "\n")

# --- 4. Build binary presence-absence matrix ----------------------------------
tf_matrix <- matrix(
  0,
  nrow     = length(accessions_dict),
  ncol     = length(all_TFs),
  dimnames = list(names(accessions_dict), all_TFs)
)

for (genome in names(accessions_dict)) {
  for (tf in accessions_dict[[genome]]) {
    tf_matrix[genome, tf] <- 1
  }
}

tf_df <- as.data.frame(tf_matrix) %>%
  rownames_to_column("Genome")

cat("Matrix:", nrow(tf_df), "genomes x", ncol(tf_df) - 1, "TF families\n")
cat("Mean TF families per genome:", round(mean(rowSums(tf_matrix)), 2), "\n")

# --- 5. Save ------------------------------------------------------------------
output_file <- file.path(output_dir, "TF_matrix.xlsx")
write_xlsx(tf_df, output_file)
cat("Saved to:", output_file, "\n")
