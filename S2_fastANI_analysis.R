#!/usr/bin/env Rscript
# =============================================================================
# S2_fastANI_analysis.R
#
# fastANI analysis: reciprocal ANI averaging, heatmaps, distribution overlay,
# descriptive statistics, ANI gap testing, clustered matrix export.
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(diptest)
  library(writexl)
})

# --- 1. Parse fastANI classical output with reciprocal averaging -----------------

read_fastani_classical <- function(file_path) {

  cat("Reading fastANI classical output with reciprocal averaging...\n")

  data <- read.table(file_path, header = FALSE, sep = "\t",
                     col.names = c("query", "reference", "ani",
                                   "coverage", "fragments"))

  data$query_clean     <- gsub("\\.fasta$", "", basename(data$query))
  data$reference_clean <- gsub("\\.fasta$", "", basename(data$reference))

  all_samples <- unique(c(data$query_clean, data$reference_clean))
  n_samples   <- length(all_samples)

  cat("Found",      n_samples,    "unique samples\n")
  cat("Processing", nrow(data),   "pairwise comparisons\n")

  # Collect reciprocal pairs
  reciprocal_data <- data.frame(
    sample1  = character(), sample2  = character(),
    ani_1to2 = numeric(),   ani_2to1 = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(data)) {
    s1      <- data$query_clean[i]
    s2      <- data$reference_clean[i]
    ani_val <- data$ani[i]
    if (s1 == s2) next

    pair   <- sort(c(s1, s2))
    std_s1 <- pair[1]; std_s2 <- pair[2]

    existing_row <- which(reciprocal_data$sample1 == std_s1 &
                          reciprocal_data$sample2 == std_s2)

    if (length(existing_row) == 0) {
      reciprocal_data <- rbind(reciprocal_data, data.frame(
        sample1  = std_s1, sample2 = std_s2,
        ani_1to2 = ifelse(s1 == std_s1, ani_val, NA),
        ani_2to1 = ifelse(s1 == std_s2, ani_val, NA),
        stringsAsFactors = FALSE
      ))
    } else {
      if (s1 == std_s1) {
        reciprocal_data$ani_1to2[existing_row] <- ani_val
      } else {
        reciprocal_data$ani_2to1[existing_row] <- ani_val
      }
    }
  }

  # Reciprocal statistics
  both     <-  !is.na(reciprocal_data$ani_1to2) &  !is.na(reciprocal_data$ani_2to1)
  one_only <- xor(!is.na(reciprocal_data$ani_1to2),
                  !is.na(reciprocal_data$ani_2to1))

  cat("\nReciprocal analysis:\n")
  cat("  Pairs with both directions:",     sum(both),               "\n")
  cat("  Pairs with one direction only:",  sum(one_only),           "\n")
  cat("  Total unique pairs:",             nrow(reciprocal_data),   "\n")

  # Average: mean where both available, single value otherwise
  reciprocal_data$average_ani       <- NA
  reciprocal_data$average_ani[both] <-
    (reciprocal_data$ani_1to2[both] + reciprocal_data$ani_2to1[both]) / 2

  one_1to2 <-  !is.na(reciprocal_data$ani_1to2) &  is.na(reciprocal_data$ani_2to1)
  one_2to1 <-   is.na(reciprocal_data$ani_1to2) & !is.na(reciprocal_data$ani_2to1)
  reciprocal_data$average_ani[one_1to2] <- reciprocal_data$ani_1to2[one_1to2]
  reciprocal_data$average_ani[one_2to1] <- reciprocal_data$ani_2to1[one_2to1]

  if (sum(both) > 0) {
    diffs <- abs(reciprocal_data$ani_1to2[both] - reciprocal_data$ani_2to1[both])
    cat("\nReciprocal difference statistics:\n")
    cat("  Mean:",   round(mean(diffs),   4), "%\n")
    cat("  Max:",    round(max(diffs),    4), "%\n")
    cat("  Median:", round(median(diffs), 4), "%\n")
  }

  # Build symmetric matrix
  ani_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples,
                       dimnames = list(all_samples, all_samples))
  diag(ani_matrix) <- 100

  for (i in 1:nrow(reciprocal_data)) {
    s1_idx <- which(all_samples == reciprocal_data$sample1[i])
    s2_idx <- which(all_samples == reciprocal_data$sample2[i])
    avg    <- reciprocal_data$average_ani[i]
    ani_matrix[s1_idx, s2_idx] <- avg
    ani_matrix[s2_idx, s1_idx] <- avg
  }

  cat("\nMatrix dimensions:", n_samples, "x", n_samples, "\n")
  return(ani_matrix)
}

# --- 2. ANI heatmap (standard resolution) ----------------------------------------

create_fastani_heatmap <- function(classical_file, output_file = NULL,
                                   width = 16, height = 16, dpi = 300,
                                   title = NULL, ani_matrix = NULL) {

  if (is.null(ani_matrix)) ani_matrix <- read_fastani_classical(classical_file)

  min_ani <- min(ani_matrix, na.rm = TRUE)
  ani_matrix[is.na(ani_matrix)] <- min_ani

  cat("ANI range:", round(min(ani_matrix), 2), "-",
                    round(max(ani_matrix), 2), "\n")

  dist_matrix <- as.dist(100 - ani_matrix)
  hc          <- hclust(dist_matrix, method = "complete")

  ani_range <- range(ani_matrix, na.rm = TRUE)
  col_fun   <- colorRamp2(
    c(ani_range[1], 98.0, 98.5, 99.0, 99.5, ani_range[2]),
    c("#000060",   "#2d0040", "#4a0080", "#DC143C", "#FFA500", "#FFFF00")
  )

  if (is.null(title)) title <- paste("FastANI Heatmap -",
                                     nrow(ani_matrix), "Strains")

  ht <- Heatmap(
    ani_matrix, name = "ANI (%)", col = col_fun,
    cluster_rows = hc, cluster_columns = hc,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp    = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7, rot = 90),
    column_title    = title,
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(
      title          = "ANI (%)",
      title_gp       = gpar(fontsize = 9, fontface = "bold"),
      labels_gp      = gpar(fontsize = 8),
      legend_height  = unit(4,   "cm"),
      legend_width   = unit(0.5, "cm")
    ),
    rect_gp           = gpar(col = "white", lwd = 0.05),
    row_dend_width    = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    row_dend_gp       = gpar(lwd = 0.8),
    column_dend_gp    = gpar(lwd = 0.8),
    width             = unit(10, "cm"),
    height            = unit(10, "cm"),
    row_names_max_width    = unit(2.5, "cm"),
    column_names_max_height = unit(2.5, "cm")
  )

  if (!is.null(output_file)) {
    png(output_file,
        width  = width  * dpi / 2.54,
        height = height * dpi / 2.54,
        res    = dpi)
    draw(ht, heatmap_legend_side = "right",
         padding = unit(c(10, 15, 10, 10), "mm"))
    dev.off()

    svg_file <- sub("\\.png$", ".svg", output_file)
    svg(svg_file, width = width / 2.54, height = height / 2.54)
    draw(ht, heatmap_legend_side = "right",
         padding = unit(c(10, 15, 10, 10), "mm"))
    dev.off()

    cat("Heatmap saved as:", output_file, "and", svg_file, "\n")
  }

  draw(ht, heatmap_legend_side = "right",
       padding = unit(c(10, 15, 10, 10), "mm"))
  return(list(matrix = ani_matrix, hclust = hc, heatmap = ht))
}

# --- 3. ANI heatmap (high resolution, for large datasets) ------------------------

create_fastani_heatmap_highres <- function(classical_file, output_file = NULL,
                                           base_width = 40, base_height = 40,
                                           dpi = 600, title = NULL,
                                           font_size = 6, ani_matrix = NULL) {

  if (is.null(ani_matrix)) ani_matrix <- read_fastani_classical(classical_file)

  min_ani <- min(ani_matrix, na.rm = TRUE)
  ani_matrix[is.na(ani_matrix)] <- min_ani
  n_samples <- nrow(ani_matrix)

  # Dynamic sizing
  sample_size_factor <- max(0.12, min(0.25, 25 / n_samples))
  dynamic_width      <- max(base_width,  n_samples * sample_size_factor)
  dynamic_height     <- max(base_height, n_samples * sample_size_factor)

  cat("Output dimensions:", round(dynamic_width, 1), "x",
                            round(dynamic_height, 1), "cm at", dpi, "DPI\n")

  dist_matrix <- as.dist(100 - ani_matrix)
  hc          <- hclust(dist_matrix, method = "complete")

  ani_range <- range(ani_matrix, na.rm = TRUE)
  col_fun   <- colorRamp2(
    c(ani_range[1], 98.0, 98.5, 99.0, 99.5, ani_range[2]),
    c("#000060",   "#2d0040", "#4a0080", "#DC143C", "#FFA500", "#FFFF00")
  )

  if (is.null(title)) title <- paste("FastANI Heatmap -",
                                     n_samples, "Strains")

  base_font   <- font_size
  title_font  <- max(14, base_font * 2)
  legend_font <- max(10, base_font * 1.5)

  ht <- Heatmap(
    ani_matrix, name = "ANI (%)", col = col_fun,
    cluster_rows = hc, cluster_columns = hc,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp    = gpar(fontsize = base_font),
    column_names_gp = gpar(fontsize = base_font, rot = 45),
    column_title    = title,
    column_title_gp = gpar(fontsize = title_font, fontface = "bold"),
    heatmap_legend_param = list(
      title         = "ANI (%)",
      title_gp      = gpar(fontsize = legend_font, fontface = "bold"),
      labels_gp     = gpar(fontsize = base_font + 2),
      legend_height = unit(8,   "cm"),
      legend_width  = unit(1.2, "cm"),
      border        = "black"
    ),
    rect_gp           = gpar(col = "white", lwd = 0.1),
    row_dend_width    = unit(max(4, dynamic_width  * 0.1), "cm"),
    column_dend_height = unit(max(4, dynamic_height * 0.1), "cm"),
    row_dend_gp       = gpar(lwd = 1.5),
    column_dend_gp    = gpar(lwd = 1.5),
    width             = unit(dynamic_width  * 0.7, "cm"),
    height            = unit(dynamic_height * 0.7, "cm"),
    row_names_max_width     = unit(max(6, dynamic_width  * 0.15), "cm"),
    column_names_max_height = unit(max(6, dynamic_height * 0.15), "cm")
  )

  if (!is.null(output_file)) {
    png(output_file,
        width  = dynamic_width  * dpi / 2.54,
        height = dynamic_height * dpi / 2.54,
        res    = dpi)
    draw(ht, heatmap_legend_side = "right",
         padding = unit(c(20, 25, 20, 20), "mm"))
    dev.off()

    file_size_mb <- round(file.info(output_file)$size / (1024^2), 1)
    cat("High-resolution heatmap saved as:", output_file,
        "(", file_size_mb, "MB)\n")
  }

  draw(ht, heatmap_legend_side = "right",
       padding = unit(c(20, 25, 20, 20), "mm"))
  return(list(matrix = ani_matrix, hclust = hc, heatmap = ht))
}

# --- 4. ANI distribution (histogram + density overlay, SVG) ----------------------

create_ani_distribution <- function(classical_file,
                                    output_file = "ani_distribution.svg",
                                    ani_matrix  = NULL) {

  if (is.null(ani_matrix)) ani_matrix <- read_fastani_classical(classical_file)

  ani_values <- ani_matrix[upper.tri(ani_matrix)]
  ani_values <- ani_values[!is.na(ani_values)]

  # Histogram uses after_stat(density) so bar heights are comparable to
  # the overlaid kernel density curve on a shared y-axis.
  p <- ggplot(data.frame(ANI = ani_values), aes(x = ANI)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 50, fill = "lightblue", alpha = 0.7,
                   color = "darkblue", linewidth = 0.3) +
    geom_density(fill = "lightcoral", alpha = 0.45,
                 color = "navy", linewidth = 1.2) +
    labs(title = paste0("ANI Distribution (",
                        length(ani_values), " comparisons)"),
         x = "ANI (%)", y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))

  # SVG output
  svg_file <- sub("\\.(png|pdf|jpg|jpeg|tif|tiff)$", ".svg",
                  output_file, ignore.case = TRUE)
  if (!grepl("\\.svg$", svg_file, ignore.case = TRUE)) {
    svg_file <- paste0(svg_file, ".svg")
  }

  ggsave(svg_file, p, width = 14, height = 9, units = "cm")
  cat("Distribution plot saved as:", svg_file, "\n")

  return(list(plot = p, values = ani_values))
}

# --- 5. Detailed statistics ------------------------------------------------------

create_detailed_statistics <- function(classical_file,
                                       output_file = "fastani_detailed_stats.txt",
                                       ani_matrix  = NULL) {

  if (is.null(ani_matrix)) ani_matrix <- read_fastani_classical(classical_file)

  ani_values <- ani_matrix[upper.tri(ani_matrix)]
  ani_values <- ani_values[!is.na(ani_values)]

  report_content <- capture.output({

    cat("FASTANI DISTRIBUTION ANALYSIS\n")
    cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Reciprocal averaging applied.\n\n")

    cat("BASIC STATISTICS:\n")
    cat("  Minimum ANI:", round(min(ani_values),    4), "%\n")
    cat("  Maximum ANI:", round(max(ani_values),    4), "%\n")
    cat("  Mean ANI:",    round(mean(ani_values),   4), "%\n")
    cat("  Median ANI:",  round(median(ani_values), 4), "%\n")
    cat("  SD:",          round(sd(ani_values),     4), "%\n")
    cat("  Q1:",          round(quantile(ani_values, 0.25), 4), "%\n")
    cat("  Q3:",          round(quantile(ani_values, 0.75), 4), "%\n")
    cat("  IQR:",         round(IQR(ani_values),    4), "%\n")
    cat("  Total pairwise comparisons:", length(ani_values), "\n")
    cat("  Strains:",                    nrow(ani_matrix),   "\n\n")

    # Peak detection
    cat("PEAK DETECTION (threshold: 3% of max density):\n")
    density_data <- density(ani_values, n = 512)
    y_vals       <- density_data$y
    x_vals       <- density_data$x
    min_height   <- max(y_vals) * 0.03

    peaks_x <- c(); peaks_y <- c()
    for (i in 2:(length(y_vals) - 1)) {
      if (y_vals[i] > y_vals[i - 1] && y_vals[i] > y_vals[i + 1] &&
          y_vals[i] > min_height) {
        peaks_x <- c(peaks_x, x_vals[i])
        peaks_y <- c(peaks_y, y_vals[i])
      }
    }

    if (length(peaks_x) > 0) {
      peak_order <- order(peaks_y, decreasing = TRUE)
      cat("  Peaks detected:", length(peaks_x), "\n")
      for (j in seq_along(peak_order)) {
        idx   <- peak_order[j]
        label <- ifelse(j == 1, "Global max", paste0("Local ", j - 1))
        cat(sprintf("  %s: %.4f%% (density: %.6f)\n",
                    label, peaks_x[idx], peaks_y[idx]))
      }
    } else {
      cat("  No peaks detected.\n")
    }

    # Distribution by similarity range
    total <- length(ani_values)
    cat("\nDISTRIBUTION BY SIMILARITY RANGE:\n")
    ranges <- list(
      ">99.5%"      = sum(ani_values >  99.5),
      "99.0-99.5%"  = sum(ani_values >  99.0 & ani_values <= 99.5),
      "98.5-99.0%"  = sum(ani_values >  98.5 & ani_values <= 99.0),
      "98.0-98.5%"  = sum(ani_values >  98.0 & ani_values <= 98.5),
      "<98.0%"      = sum(ani_values <= 98.0)
    )
    for (nm in names(ranges)) {
      cat(sprintf("  %s: %d (%.2f%%)\n",
                  nm, ranges[[nm]], ranges[[nm]] / total * 100))
    }

    # ANI gap ranges
    total_range <- max(ani_values) - min(ani_values)
    cat("\nANI GAP RANGES (observed/expected ratio):\n")
    gap_defs <- list(
      "99.0-99.5%" = c(99.0, 99.5),
      "99.2-99.8%" = c(99.2, 99.8),
      "99.5-99.8%" = c(99.5, 99.8),
      "99.6-99.8%" = c(99.6, 99.8)
    )
    for (nm in names(gap_defs)) {
      lo    <- gap_defs[[nm]][1]
      hi    <- gap_defs[[nm]][2]
      n_in  <- sum(ani_values >= lo & ani_values <= hi)
      exp_n <- round(total * (hi - lo) / total_range)
      cat(sprintf("  %s: %d obs / %d exp (ratio: %.3f)\n",
                  nm, n_in, exp_n,
                  round(n_in / max(exp_n, 1), 3)))
    }

    # Genomovar analysis
    cat("\nGENOMOVAR ANALYSIS (ANI > 99.5%):\n")
    genomovar_n <- sum(ani_values > 99.5)
    cat(sprintf("  Same genomovar (>99.5%%): %d (%.2f%%)\n",
                genomovar_n, genomovar_n / total * 100))

  }, type = "output")

  writeLines(report_content, output_file)
  cat("Statistics saved to:", output_file, "\n")
  cat(paste(report_content, collapse = "\n"))
  return(report_content)
}

# --- 6. Statistical testing of ANI gap -------------------------------------------

test_ani_gap <- function(classical_file,
                         output_file = "ani_gap_statistics.txt",
                         gap_ranges  = list(
                           "99.0-99.5%" = c(99.0, 99.5),
                           "99.2-99.8%" = c(99.2, 99.8),
                           "99.5-99.8%" = c(99.5, 99.8),
                           "99.6-99.8%" = c(99.6, 99.8)
                         ),
                         n_bootstrap = 1000, ani_matrix = NULL) {

  if (is.null(ani_matrix)) ani_matrix <- read_fastani_classical(classical_file)

  ani_values  <- ani_matrix[upper.tri(ani_matrix)]
  ani_values  <- ani_values[!is.na(ani_values)]
  total       <- length(ani_values)
  total_range <- max(ani_values) - min(ani_values)

  dip_result <- NULL   # captured for return value

  report_content <- capture.output({

    cat("ANI GAP STATISTICAL ANALYSIS\n")
    cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

    # --- Hartigan's dip test ---
    cat("HARTIGAN'S DIP TEST FOR MULTIMODALITY:\n")
    cat("  H0: unimodal distribution\n")
    cat("  H1: multimodal distribution\n\n")

    dip_result <<- dip.test(ani_values)
    cat("  Dip statistic (D):", round(dip_result$statistic, 6), "\n")
    cat("  P-value:",           format(dip_result$p.value,
                                       scientific = TRUE, digits = 4), "\n")

    if (dip_result$p.value < 0.001) {
      cat("  Interpretation: Highly significant (p < 0.001).\n")
    } else if (dip_result$p.value < 0.05) {
      cat("  Interpretation: Significant (p < 0.05).\n")
    } else {
      cat("  Interpretation: Not significant - cannot reject unimodality.\n")
    }

    # --- Gap underrepresentation for each range ---
    for (range_name in names(gap_ranges)) {
      gap_low  <- gap_ranges[[range_name]][1]
      gap_high <- gap_ranges[[range_name]][2]

      cat(sprintf("\nANI GAP ANALYSIS (%s):\n", range_name))

      below  <- sum(ani_values < gap_low)
      within <- sum(ani_values >= gap_low & ani_values <= gap_high)
      above  <- sum(ani_values > gap_high)

      gap_width  <- gap_high - gap_low
      exp_prop   <- gap_width / total_range
      exp_n      <- round(total * exp_prop)

      cat(sprintf("  Below: %d | Within: %d | Above: %d\n",
                  below, within, above))
      cat(sprintf("  Expected if uniform: %d | Obs/Exp ratio: %.3f\n",
                  exp_n, within / max(exp_n, 1)))

      # Chi-square test
      observed   <- c(below + above, within)
      expected_p <- c(1 - exp_prop, exp_prop)

      if (all(total * expected_p > 5)) {
        chi_result <- chisq.test(observed, p = expected_p)
        cat(sprintf("  Chi-square: %.2f, p = %s\n",
                    chi_result$statistic,
                    format(chi_result$p.value,
                           scientific = TRUE, digits = 4)))
        if (chi_result$p.value < 0.05 && within < exp_n) {
          cat("  -> Underrepresentation is statistically significant.\n")
        }
      } else {
        cat("  Note: Expected counts too low for chi-square test.\n")
      }
    }

    # --- Density-based gap detection ---
    cat("\nDENSITY-BASED GAP DETECTION (99.0-99.9% range):\n")
    density_data <- density(ani_values, n = 1024)
    high_idx     <- which(density_data$x >= 99.0 & density_data$x <= 99.9)

    if (length(high_idx) > 0) {
      min_idx     <- high_idx[which.min(density_data$y[high_idx])]
      gap_center  <- density_data$x[min_idx]
      gap_density <- density_data$y[min_idx]
      max_density <- max(density_data$y)

      cat(sprintf("  Density minimum: %.3f%% ANI (%.2f%% of max density)\n",
                  gap_center, gap_density / max_density * 100))
    }

    # --- Bootstrap peak stability ---
    cat(sprintf("\nBOOTSTRAP PEAK STABILITY (B = %d):\n", n_bootstrap))

    # Auto-detect peaks from data
    d            <- density(ani_values, n = 512)
    pk_x         <- c()
    pk_y         <- c()
    pk_threshold <- max(d$y) * 0.03
    for (i in 2:(length(d$y) - 1)) {
      if (d$y[i] > d$y[i - 1] && d$y[i] > d$y[i + 1] &&
          d$y[i] > pk_threshold) {
        pk_x <- c(pk_x, d$x[i])
        pk_y <- c(pk_y, d$y[i])
      }
    }
    peak_targets <- round(pk_x, 2)
    cat("  Peak targets (auto-detected):",
        paste(peak_targets, collapse = ", "), "\n")

    find_peak_near <- function(x, target, tol = 0.02,
                               min_height_frac = 0.03) {
      d     <- density(x)
      peaks <- which(diff(sign(diff(d$y))) == -2) + 1
      peaks <- peaks[d$y[peaks] > max(d$y) * min_height_frac]
      any(abs(d$x[peaks] - target) < tol)
    }

    set.seed(1)
    for (pt in peak_targets) {
      freq <- mean(replicate(n_bootstrap, {
        find_peak_near(sample(ani_values, replace = TRUE), pt)
      }))
      cat(sprintf("  Peak at %.2f%%: recovery = %.3f\n", pt, freq))
    }

    # --- Genomovar ---
    cat("\nGENOMOVAR ANALYSIS (ANI > 99.5%):\n")
    genomovar_n <- sum(ani_values > 99.5)
    cat(sprintf("  Same genomovar: %d (%.2f%%)\n",
                genomovar_n, genomovar_n / total * 100))

  }, type = "output")

  writeLines(report_content, output_file)
  cat("Gap analysis saved to:", output_file, "\n")
  cat(paste(report_content, collapse = "\n"))

  return(list(dip_test = dip_result, total = total))
}

# --- 7. Export clustered ANI matrix to Excel -------------------------------------

export_for_manual_review <- function(classical_file,
                                     output_file = "clustered_matrix.xlsx",
                                     ani_matrix  = NULL) {

  if (is.null(ani_matrix)) ani_matrix <- read_fastani_classical(classical_file)

  min_ani <- min(ani_matrix, na.rm = TRUE)
  ani_matrix[is.na(ani_matrix)] <- min_ani

  dist_matrix <- as.dist(100 - ani_matrix)
  hc          <- hclust(dist_matrix, method = "complete")

  clustered_order  <- hc$order
  clustered_matrix <- ani_matrix[clustered_order, clustered_order]
  strain_names     <- rownames(ani_matrix)[clustered_order]

  clustered_df              <- as.data.frame(round(clustered_matrix, 3))
  rownames(clustered_df)    <- strain_names
  colnames(clustered_df)    <- strain_names

  clustered_df_with_names <- data.frame(
    Strain = strain_names, clustered_df,
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # Cluster assignments at multiple cut levels
  cut_levels   <- c(3, 5, 8, 10, 15, 20)
  cluster_cuts <- list()

  for (k in cut_levels) {
    clusters <- cutree(hc, k = k)
    cluster_cuts[[paste0("Clusters_", k)]] <- data.frame(
      Position = 1:length(strain_names),
      Strain   = strain_names,
      Cluster  = clusters[clustered_order],
      stringsAsFactors = FALSE
    )
  }

  # Cluster boundary guide
  cluster_guide <- data.frame(
    Cut_Level    = integer(),  Cluster     = integer(),
    Start        = integer(),  End         = integer(),
    Count        = integer(),  First_Strain = character(),
    Last_Strain  = character(),
    stringsAsFactors = FALSE
  )

  for (k in cut_levels) {
    ordered_clusters <- cutree(hc, k = k)[clustered_order]
    for (cl in 1:k) {
      pos <- which(ordered_clusters == cl)
      if (length(pos) > 0) {
        cluster_guide <- rbind(cluster_guide, data.frame(
          Cut_Level    = k,                    Cluster = cl,
          Start        = min(pos),             End     = max(pos),
          Count        = length(pos),
          First_Strain = strain_names[min(pos)],
          Last_Strain  = strain_names[max(pos)],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  excel_sheets <- c(
    list(ANI_Matrix = clustered_df_with_names,
         Cluster_Guide = cluster_guide),
    cluster_cuts
  )

  write_xlsx(excel_sheets, output_file)
  cat("Clustered matrix exported to:", output_file, "\n")
  cat("  Matrix:",  nrow(clustered_matrix), "x",
                   ncol(clustered_matrix), "\n")
  cat("  Sheets: ANI_Matrix, Cluster_Guide, Clusters_3/5/8/10/15/20\n")

  return(list(matrix = clustered_matrix,
              strain_names = strain_names,
              hclust = hc))
}

# --- 8. Complete analysis wrapper ------------------------------------------------


create_complete_fastani_analysis <- function(classical_file,
                                             output_prefix = "fastani_analysis",
                                             highres       = FALSE,
                                             export_excel  = FALSE) {

  cat("Running complete FastANI analysis...\n\n")

  # Parse once, reuse
  ani_matrix <- read_fastani_classical(classical_file)

  # Heatmap
  if (highres) {
    heatmap_result <- create_fastani_heatmap_highres(
      classical_file,
      output_file = paste0(output_prefix, "_heatmap_highres.png"),
      ani_matrix  = ani_matrix
    )
  } else {
    heatmap_result <- create_fastani_heatmap(
      classical_file,
      output_file = paste0(output_prefix, "_heatmap.png"),
      ani_matrix  = ani_matrix
    )
  }

  # Distribution (single overlay plot, SVG)
  dist_result <- create_ani_distribution(
    classical_file,
    output_file = paste0(output_prefix, "_distribution.svg"),
    ani_matrix  = ani_matrix
  )

  # Statistics
  stats_result <- create_detailed_statistics(
    classical_file,
    output_file = paste0(output_prefix, "_detailed_stats.txt"),
    ani_matrix  = ani_matrix
  )

  # Gap analysis
  gap_result <- test_ani_gap(
    classical_file,
    output_file = paste0(output_prefix, "_gap_statistics.txt"),
    ani_matrix  = ani_matrix
  )

  # Optional Excel export
  excel_result <- NULL
  if (export_excel) {
    excel_result <- export_for_manual_review(
      classical_file,
      output_file = paste0(output_prefix, "_clustered_matrix.xlsx"),
      ani_matrix  = ani_matrix
    )
  }

  cat("\nAnalysis complete. Output prefix:", output_prefix, "\n")
  return(list(heatmap      = heatmap_result,
              distribution = dist_result,
              statistics   = stats_result,
              gap          = gap_result,
              excel        = excel_result))
}

# --- 9. Main entry point: runs when script is executed with Rscript --------------

args <- commandArgs(trailingOnly = TRUE)

# --- EDIT DEFAULTS OR PASS AS COMMAND-LINE ARGS ----
classical_file <- if (length(args) >= 1) args[1] else "input/fastani_classical.tsv"
output_dir     <- if (length(args) >= 2) args[2] else "output"
# --------------------------------------------------

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Run standard-resolution analysis with Excel export.
# For large datasets (e.g. 313 strains), highres = TRUE.
result <- create_complete_fastani_analysis(
  classical_file,
  output_prefix = file.path(output_dir, "fastani_complete"),
  highres       = FALSE,
  export_excel  = TRUE
)

# Find maximum non-self ANI pair
ani_mat       <- result$heatmap$matrix
diag(ani_mat) <- NA
max_ani       <- max(ani_mat, na.rm = TRUE)
cat("\nMaximum non-self ANI:", max_ani, "%\n")
max_pos <- which(ani_mat == max_ani, arr.ind = TRUE)
for (i in 1:nrow(max_pos)) {
  cat("Pair:", rownames(ani_mat)[max_pos[i, 1]], "<->",
               colnames(ani_mat)[max_pos[i, 2]], "\n")
}
