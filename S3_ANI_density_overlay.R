#!/usr/bin/env Rscript
# =============================================================================
# S3_ANI_density_overlay.R
#
# Overlay kernel density plot of pairwise ANI distributions from two datasets
# (40-strain local panel vs 313-strain species-wide panel).
#
# Author: Eva Stare
# =============================================================================


# --- 1. Libraries -------------------------------------------------------------
library(tidyverse)
library(svglite)


# --- 2. User configuration ----------------------------------------------------
# Paths to FastANI classical-format output files.
# Each file should contain reciprocal pairs (both A-->B and B-->A), which are
# averaged below to produce symmetric per-pair ANI values.
classical_file_40  <- "input/fastani_40_strain_output.txt"
classical_file_313 <- "input/fastani_313_strain_output.txt"

# Directory for figure output (created if missing)
output_dir <- "figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Colors
color_40  <- "#0072B2"   # blue
color_313 <- "#D55E00"   # vermillion


# --- 3. Helper functions ------------------------------------------------------

#' Read FastANI classical-format output and build a symmetric ANI matrix
#'
#' Reciprocal ANI values (A->B and B->A) are averaged per strain pair.
#' Self-comparisons are set to 100.
#'
#' @param file_path Path to FastANI classical-format tab-separated file.
#' @return Symmetric numeric matrix of ANI values (rows/cols = strain names).
read_fastani_classical <- function(file_path) {

  data <- read.table(file_path, header = FALSE, sep = "\t",
                     stringsAsFactors = FALSE)
  colnames(data)[1:3] <- c("query", "reference", "ANI")

  # Strip paths and file extensions from strain names
  data$query     <- tools::file_path_sans_ext(basename(data$query))
  data$reference <- tools::file_path_sans_ext(basename(data$reference))

  all_samples <- sort(unique(c(data$query, data$reference)))
  n_samples   <- length(all_samples)

  m <- matrix(NA_real_, nrow = n_samples, ncol = n_samples,
              dimnames = list(all_samples, all_samples))
  diag(m) <- 100

  for (i in seq_len(nrow(data))) {
    s1 <- data$query[i]
    s2 <- data$reference[i]
    if (s1 == s2) next
    existing <- m[s1, s2]
    new_val  <- data$ANI[i]
    m[s1, s2] <- if (is.na(existing)) new_val else mean(c(existing, new_val))
    m[s2, s1] <- m[s1, s2]
  }

  cat(sprintf("  Loaded %d x %d ANI matrix from: %s\n",
              n_samples, n_samples, basename(file_path)))
  m
}


#' Extract unique pairwise ANI values from a symmetric matrix (upper triangle)
get_pairwise_values <- function(ani_matrix) {
  vals <- ani_matrix[upper.tri(ani_matrix)]
  vals[!is.na(vals)]
}


# --- 4. Load and summarise data -----------------------------------------------
message("Reading FastANI output files...")
matrix_40      <- read_fastani_classical(classical_file_40)
ani_values_40  <- get_pairwise_values(matrix_40)

matrix_313     <- read_fastani_classical(classical_file_313)
ani_values_313 <- get_pairwise_values(matrix_313)

cat(sprintf("\n40-strain panel:  %s pairs (ANI range %.2f - %.2f%%)\n",
            format(length(ani_values_40), big.mark = ","),
            min(ani_values_40), max(ani_values_40)))
cat(sprintf("313-strain panel: %s pairs (ANI range %.2f - %.2f%%)\n\n",
            format(length(ani_values_313), big.mark = ","),
            min(ani_values_313), max(ani_values_313)))


# --- 5. Build combined data frame ---------------------------------------------
label_40  <- sprintf("40-strain panel (n = %s)",
                     format(length(ani_values_40),  big.mark = ","))
label_313 <- sprintf("313-strain panel (n = %s)",
                     format(length(ani_values_313), big.mark = ","))

ani_df <- bind_rows(
  tibble(ANI = ani_values_40,  dataset = label_40),
  tibble(ANI = ani_values_313, dataset = label_313)
) %>%
  mutate(dataset = factor(dataset, levels = c(label_40, label_313)))


# --- 6. Figure 2C: overlay density plot ---------------------------------------
p_overlay <- ggplot(ani_df, aes(x = ANI, fill = dataset, color = dataset)) +
  # Shade the 99.0–99.5% ANI gap region
  annotate("rect", xmin = 99.0, xmax = 99.5,
           ymin = 0, ymax = Inf,
           alpha = 0.15, fill = "grey40") +
  geom_density(alpha = 0.4, linewidth = 0.9) +
  scale_x_continuous(
    breaks       = seq(97, 100, by = 0.5),
    minor_breaks = seq(97, 100, by = 0.1),
    limits       = c(97, 100),
    expand       = c(0, 0)
  ) +
  scale_fill_manual(values  = c(color_40, color_313)) +
  scale_color_manual(values = c(color_40, color_313)) +
  labs(x = "ANI (%)", y = "Density", fill = NULL, color = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title       = element_text(size = 16, face = "bold"),
    axis.text        = element_text(size = 12),
    legend.text      = element_text(size = 12),
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )


# --- 7. Save figure -----------------------------------------------------------
ggsave(file.path(output_dir, "Figure2C_ANI_density_overlay.svg"),
       plot = p_overlay, width = 8, height = 5, device = svglite)
ggsave(file.path(output_dir, "Figure2C_ANI_density_overlay.png"),
       plot = p_overlay, width = 8, height = 5, dpi = 300)

cat("Figure saved to: ", output_dir, "\n", sep = "")
