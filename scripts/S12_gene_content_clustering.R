#!/usr/bin/env Rscript
# =============================================================================
# S12_gene_content_clustering.R
#
# Accessory gene content clustering via Jaccard distance + UPGMA.
#
# Constructs a UPGMA dendrogram from binary gene presence/absence data using
# Jaccard distance, for visualisation of strain relationships based on
# accessory gene content.
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(vegan)
  library(ape)
})

# --- 1. Configuration ---------------------------------------------------------
feature_file <- "input/Mega_matrix.xlsx"
sheet_name   <- "Filtered_Genes"           # adjust to your sheet name
output_dir   <- "output"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. Load and prepare data -------------------------------------------------
df <- as.data.frame(read_excel(feature_file, sheet = sheet_name),
                    check.names = FALSE)
strain_names     <- as.character(df[[1]])
mat              <- df[, -1, drop = FALSE]
rownames(mat)    <- strain_names

# Convert to integer and validate
mat[] <- lapply(mat, function(x) suppressWarnings(as.integer(x)))
if (anyNA(mat)) stop("NAs introduced during conversion")

cat("Loaded:", nrow(mat), "strains x", ncol(mat), "genes\n")

# --- 3. Remove invariant genes ------------------------------------------------
var_cols <- vapply(mat, function(x) length(unique(x)) > 1, logical(1))
mat      <- mat[, var_cols, drop = FALSE]

cat("Variable genes retained:", ncol(mat),       "\n")
cat("Invariant genes removed:", sum(!var_cols),  "\n")

# --- 4. Compute Jaccard distances and cluster ---------------------------------
mat_num              <- as.matrix(mat)
storage.mode(mat_num) <- "numeric"

dist_jaccard <- vegdist(mat_num, method = "jaccard", binary = TRUE)
hc           <- hclust(dist_jaccard, method = "average")

cat("Clustering complete:", length(hc$labels), "strains\n")

# --- 5. Export Newick tree ----------------------------------------------------
phy <- as.phylo(hc)
phy$edge.length <- phy$edge.length * 2   # correct for as.phylo() halving

tree_file <- file.path(output_dir, "gene_content_UPGMA_rescaled.nwk")
write.tree(phy, file = tree_file)

cat("Newick tree saved to:", tree_file, "\n")
cat("Branch lengths represent Jaccard distances.\n")
