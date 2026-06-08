#!/usr/bin/env Rscript
# =============================================================================
# S4_pangenome_plots.R
#
# Pangenome accumulation and BLAST identity plots from Roary output.
# Adapted from create_pan_genome_plots.R (Page et al., 2015; Roary v3.13.0;
# https://github.com/sanger-pathogens/Roary).

# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# --- 1. Arguments & directories -----------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "input"
output_dir <- if (length(args) >= 2) args[2] else "output_plots"

if (!dir.exists(input_dir)) {
  stop("Input directory '", input_dir, "' does not exist.")
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. Locate Roary output files ---------------------------------------------
find_one <- function(pattern) {
  hits <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) stop("No file matching '", pattern, "' in ", input_dir)
  if (length(hits) > 1)  stop("Multiple files matching '", pattern, "' in ", input_dir)
  hits
}

files <- list(
  new       = find_one(".*Number of New Genes.*\\.Rtab$"),
  conserved = find_one(".*Number of Conserved Genes.*\\.Rtab$"),
  pan       = find_one(".*Number of Genes in Pan Geneome.*\\.Rtab$"),
  unique    = find_one(".*Number of Unique Genes.*\\.Rtab$"),
  blast     = find_one(".*Blast Identity Frequencies.*\\.Rtab$")
)

# --- 3. Boxplots: accumulation curves -----------------------------------------
save_boxplot <- function(path, data, title) {
  pdf(file.path(output_dir, path))
  boxplot(data,
          main      = title,
          xlab      = "No. of genomes",
          ylab      = "No. of genes",
          varwidth  = TRUE,
          ylim      = c(0, max(data)),
          outline   = FALSE)
  dev.off()
}

data_new  <- read.table(files$new)
data_con  <- read.table(files$conserved)
data_pan  <- read.table(files$pan)
data_uniq <- read.table(files$unique)

save_boxplot("Number_of_new_genes.pdf",        data_new,  "Number of new genes")
save_boxplot("Number_of_conserved_genes.pdf",  data_con,  "Number of conserved genes")
save_boxplot("Number_of_genes_in_pan_genome.pdf", data_pan, "No. of genes in the pan-genome")
save_boxplot("Number_of_unique_genes.pdf",     data_uniq, "Number of unique genes")

# --- 4. BLAST identity frequency plot -----------------------------------------
pdf(file.path(output_dir, "Blast_identity_frequencies.pdf"))
plot(read.table(files$blast),
     main = "BLASTp hits by percentage identity",
     xlab = "BLAST percentage identity",
     ylab = "No. of BLAST hits")
dev.off()

# --- 5. Helper: mean-curve line plot saved as PDF + PNG -----------------------
plot_mean_curves <- function(series_a, series_b, label_a, label_b,
                             legend_corner, outfile_stem) {
  df <- data.frame(
    value   = c(series_a, series_b),
    genomes = rep(seq_along(series_a), 2),
    Key     = factor(rep(c(label_a, label_b), each = length(series_a)),
                     levels = c(label_a, label_b))
  )
  p <- ggplot(df, aes(x = genomes, y = value, linetype = Key)) +
    geom_line() +
    ylim(1, max(df$value)) +
    xlim(1, length(series_a)) +
    xlab("No. of genomes") +
    ylab("No. of genes") +
    theme_bw(base_size = 16) +
    theme(legend.justification = legend_corner,
          legend.position      = legend_corner)
  ggsave(file.path(output_dir, paste0(outfile_stem, ".png")), p)
  ggsave(file.path(output_dir, paste0(outfile_stem, ".pdf")), p)
}

plot_mean_curves(colMeans(data_con), colMeans(data_pan),
                 "Conserved genes", "Total genes",
                 legend_corner = c(0, 1),
                 outfile_stem  = "conserved_vs_total_genes")

plot_mean_curves(colMeans(data_uniq), colMeans(data_new),
                 "Unique genes", "New genes",
                 legend_corner = c(1, 1),
                 outfile_stem  = "unique_vs_new_genes")

message("Plots written to: ", normalizePath(output_dir))
