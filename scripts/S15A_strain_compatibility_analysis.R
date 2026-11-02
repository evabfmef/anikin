#!/usr/bin/env Rscript
# =============================================================================
# S15A_strain_compatibility_analysis.R
#
# Kin compatibility: fastANI thresholds and within-group concordance.
#
# Analyses the concordance between pairwise swarming phenotypes and fastANI
# genomic similarity across the strain panel, produces the fastANI threshold
# summary reported in the main text, generates the fastANI distribution plot,
# and stratifies nonconforming pairs into strong / moderate / borderline
# classes.
#
# CONFIGURATION
#   All paths, thresholds, phenotype definitions and plot settings are read
#   from a YAML configuration file (default: config_S15.yml) and are not
#   hard-coded below. Override individual values on the command line.
#
# USAGE
#   Rscript S15A_strain_compatibility_analysis.R [options]
#
#     --config=FILE     configuration file      (default: config_S15.yml)
#     --input=FILE      input .xlsx             (overrides config)
#     --output=DIR      output directory        (overrides config)
#     --kin-cutoff=X    phenotype >= X is kin   (overrides config)
#     --help            show this message
#
#   Examples:
#     Rscript S15A_strain_compatibility_analysis.R
#     Rscript S15A_strain_compatibility_analysis.R --config=my_config.yml
#     Rscript S15A_strain_compatibility_analysis.R \
#         --input=data/panel_B.xlsx --output=results/panel_B
#
#
# Inputs : as set in config (default data/strain_interactions.xlsx)
# Outputs: <output_dir>/within_group_summary.csv
#          <output_dir>/strong_nonconformers.csv
#          <output_dir>/moderate_nonconformers.csv
#          <output_dir>/borderline_pairs.csv
#          <output_dir>/fastani_vs_phenotype.png
#          <output_dir>/threshold_plot.png
#          <output_dir>/within_vs_between_groups.png
#          <output_dir>/fastani_distribution.png
#          <output_dir>/summary.txt
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(knitr)
  library(gridExtra)
  library(yaml)
})

# =============================================================================
# CONFIGURATION LAYER
# No data processing below this point until section 1.
# =============================================================================

usage <- function() {
  cat("
Usage: Rscript S15A_strain_compatibility_analysis.R [options]

  --config=FILE     configuration file    (default: config_S15.yml)
  --input=FILE      input .xlsx           (overrides config)
  --output=DIR      output directory      (overrides config)
  --kin-cutoff=X    phenotype >= X is kin (overrides config)
  --help            show this message

")
  quit(status = 0)
}

parse_args <- function(args) {
  if ("--help" %in% args || "-h" %in% args) usage()
  out <- list()
  for (a in args) {
    if (!grepl("^--[a-z-]+=", a)) {
      stop("Unrecognised argument: ", a,
           "\n  Run with --help for usage.", call. = FALSE)
    }
    key <- sub("^--([a-z-]+)=.*$", "\\1", a)
    val <- sub("^--[a-z-]+=", "", a)
    out[[key]] <- val
  }
  out
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

config_file <- if (!is.null(opts$config)) opts$config else "config_S15.yml"
if (!file.exists(config_file)) {
  stop("Configuration file not found: ", config_file,
       "\n  Expected config_S15.yml in the working directory, or pass",
       " --config=FILE.", call. = FALSE)
}
cfg <- yaml::read_yaml(config_file)

# Command-line arguments override configuration file values.
input_file  <- if (!is.null(opts$input))  opts$input  else cfg$paths$input_file
output_dir  <- if (!is.null(opts$output)) opts$output else cfg$paths$output_dir
kin_cutoff  <- if (!is.null(opts[["kin-cutoff"]]))
                 as.numeric(opts[["kin-cutoff"]]) else cfg$phenotype$kin_cutoff

weak_kin_cutoff <- cfg$phenotype$weak_kin_cutoff
lysis_value     <- cfg$phenotype$lysis_value

# ANI bins (descriptive; see note in header).
ani_strong     <- cfg$ani$strong_nonconformer
ani_moderate   <- cfg$ani$moderate_nonconformer
ani_border_lo  <- cfg$ani$borderline_lower
ani_border_hi  <- cfg$ani$borderline_upper

# Phenotype scale: values, labels and colours are defined once in the config
# and expanded here into the named vectors ggplot expects.
phenotype_values <- vapply(cfg$phenotype$levels,
                           function(x) as.numeric(x$value), numeric(1))
phenotype_labels <- setNames(
  vapply(cfg$phenotype$levels, function(x) x$label, character(1)),
  as.character(phenotype_values))
phenotype_colors <- setNames(
  vapply(cfg$phenotype$levels, function(x) x$color, character(1)),
  as.character(phenotype_values))
# Ordered high-to-low for consistent legend ordering.
phenotype_order <- phenotype_labels[order(phenotype_values, decreasing = TRUE)]

# Plot settings.
plot_dpi    <- cfg$plots$dpi
plot_base   <- cfg$plots$base_size
plot_bins   <- cfg$plots$histogram_bins

# Category labels are BUILT from the threshold values rather than written out,
# so that changing a threshold in the config updates every label, caption and
# summary line automatically.
lab_strong   <- sprintf("Strong nonconformers (>=%.1f%% + strong boundary)",   ani_strong)
lab_moderate <- sprintf("Moderate nonconformers (>=%.1f%% + any non-kin)",
                        ani_moderate)
lab_border   <- sprintf("Borderline pairs (%.1f-%.1f%% ANI)",
                        ani_border_lo, ani_border_hi)

if (!file.exists(input_file))
  stop("Input file not found: ", input_file,
       "\n  Set paths:input_file in ", config_file,
       " or pass --input=FILE.", call. = FALSE)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Configuration:      ", config_file, "\n", sep = "")
cat("Input:              ", input_file,  "\n", sep = "")
cat("Output directory:   ", output_dir,  "\n", sep = "")
cat(sprintf("Kin cutoff:          phenotype >= %.2f\n", kin_cutoff))

# =============================================================================
# ANALYSIS
# =============================================================================

# --- 1. Load data -------------------------------------------------------------
data <- read_excel(input_file) |>
  mutate(
    Phenotype_interaction = as.numeric(Phenotype_interaction),
    fastANI               = as.numeric(fastANI)
  )

required_cols <- c("S1", "S2", "Phenotype_interaction",
                   "Group_S1", "Group_S2", "fastANI")
missing_cols  <- setdiff(required_cols, names(data))
if (length(missing_cols))
  stop("Input file is missing required column(s): ",
       paste(missing_cols, collapse = ", "), call. = FALSE)

cat("Pairs loaded:",    nrow(data),                            "\n")
cat("Unique strains:",  length(unique(c(data$S1, data$S2))),   "\n")

# --- 2. Within-group concordance ----------------------------------------------
# For each group with >=2 tested members, summarise within-group swarming
# phenotypes and mean fastANI.
#   kin_consistent      : all phenotypes >= kin_cutoff
#   weak_kin_consistent : all phenotypes >= weak_kin_cutoff
within_group <- data |>
  filter(Group_S1 == Group_S2) |>
  group_by(group = Group_S1) |>
  summarise(
    n_pairs             = n(),
    min_phenotype       = min( Phenotype_interaction, na.rm = TRUE),
    mean_phenotype      = mean(Phenotype_interaction, na.rm = TRUE),
    mean_fastani        = mean(fastANI,               na.rm = TRUE),
    kin_consistent      = all(Phenotype_interaction >= kin_cutoff,
                              na.rm = TRUE),
    weak_kin_consistent = all(Phenotype_interaction >= weak_kin_cutoff,
                              na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(mean_phenotype))

write.csv(within_group,
          file.path(output_dir, "within_group_summary.csv"),
          row.names = FALSE)

print(kable(within_group, digits = 2,
            caption = "Within-group swarming phenotype and fastANI summary."))

# --- 3. fastANI ranges per phenotype class ------------------------------------
# The minimum observed fastANI per class is a DESCRIPTIVE boundary statistic
# of this panel: it is the most sample-dependent statistic available and is
# not a fitted decision threshold. S15B quantifies its uncertainty.
threshold_summary <- function(d, mask, label) {
  vals <- d$fastANI[mask & !is.na(d$fastANI)]
  if (!length(vals)) return(NULL)
  tibble(class = label, n = length(vals),
         min_fastani = min(vals), max_fastani = max(vals))
}

# Per-level classes are generated from the configured phenotype scale.
thr_levels <- bind_rows(lapply(seq_along(phenotype_values), function(i) {
  v <- phenotype_values[i]
  threshold_summary(data, data$Phenotype_interaction == v,
                    phenotype_labels[[as.character(v)]])
}))

thr <- bind_rows(
  thr_levels,
  threshold_summary(data, data$Phenotype_interaction >= kin_cutoff,
                    sprintf("Any kin (>= %.2g)", kin_cutoff))
)

print(kable(thr, digits = 2,
            caption = "fastANI ranges by phenotype class."))

# --- 4. Pearson correlation ---------------------------------------------------
cor_res <- cor.test(data$fastANI, data$Phenotype_interaction,
                    method = "pearson")

cor_tbl <- tibble(
  r       = cor_res$estimate,
  p_value = cor_res$p.value,
  ci_lo   = cor_res$conf.int[1],
  ci_hi   = cor_res$conf.int[2]
)
print(kable(cor_tbl, digits = c(4, 3, 4, 4),
            caption = "Pearson correlation between fastANI and phenotype."))

# --- 5. Nonconforming pairs stratification ------------------------------------
# Bins are defined in the configuration file; labels are derived from those
# values so that thresholds and text can never diverge.
strong_nonconformers <- data |>
  filter(fastANI >= ani_strong, Phenotype_interaction == lysis_value) |>
  arrange(desc(fastANI))

moderate_nonconformers <- data |>
  filter(fastANI >= ani_moderate,
         Phenotype_interaction <= weak_kin_cutoff) |>
  arrange(desc(fastANI), Phenotype_interaction)

borderline_pairs <- data |>
  filter(fastANI >= ani_border_lo, fastANI < ani_border_hi) |>
  arrange(desc(fastANI))

nonconf_counts <- tibble(
  Category = c(lab_strong, lab_moderate, lab_border),
  n = c(nrow(strong_nonconformers),
        nrow(moderate_nonconformers),
        nrow(borderline_pairs))
)
print(kable(nonconf_counts,
            caption = "Nonconforming-pair stratification (counts)."))

write.csv(strong_nonconformers,
          file.path(output_dir, "strong_nonconformers.csv"),
          row.names = FALSE)
write.csv(moderate_nonconformers,
          file.path(output_dir, "moderate_nonconformers.csv"),
          row.names = FALSE)
write.csv(borderline_pairs,
          file.path(output_dir, "borderline_pairs.csv"),
          row.names = FALSE)

if (nrow(strong_nonconformers)) {
  print(
    strong_nonconformers |>
      select(S1, S2, Phenotype_interaction, Group_S1, Group_S2, fastANI) |>
      kable(digits = 3, caption = paste("Strong nonconforming pairs:",
                                        lab_strong))
  )
}

# --- 6. Plots -----------------------------------------------------------------

# --- 6a. fastANI vs. phenotype ------------------------------------------------
p_scatter <- ggplot(data, aes(fastANI, Phenotype_interaction)) +
  geom_jitter(aes(color = factor(Phenotype_interaction)),
              size = 2.5, alpha = 0.7, height = 0.02, width = 0) +
  geom_smooth(method = "loess", se = TRUE, color = "black", alpha = 0.2) +
  scale_color_manual(values = phenotype_colors,
                     labels = phenotype_labels,
                     name   = "Phenotype") +
  labs(x = "fastANI (%)", y = "Phenotype score") +
  theme_minimal(base_size = plot_base) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fastani_vs_phenotype.png"), p_scatter,
       width  = cfg$plots$scatter$width,
       height = cfg$plots$scatter$height, dpi = plot_dpi)

# --- 6b. Threshold plot -------------------------------------------------------
# Reference lines mark the minimum OBSERVED fastANI per class. These are
# descriptive, not fitted; see S15B.
top_label      <- phenotype_labels[[as.character(max(phenotype_values))]]
strong_kin_min <- thr$min_fastani[thr$class == top_label]
any_kin_min    <- thr$min_fastani[thr$class ==
                                  sprintf("Any kin (>= %.2g)", kin_cutoff)]

p_thresh <- ggplot(data, aes(fastANI, Phenotype_interaction)) +
  geom_point(aes(color = factor(Phenotype_interaction)),
             size = 2.5, alpha = 0.7) +
  { if (length(strong_kin_min))
      geom_vline(xintercept = strong_kin_min, linetype = "dashed",
                 color = "red", linewidth = 1) } +
  { if (length(any_kin_min))
      geom_vline(xintercept = any_kin_min, linetype = "dashed",
                 color = "darkgreen", linewidth = 1) } +
  scale_color_manual(values = phenotype_colors,
                     labels = phenotype_labels,
                     name   = "Phenotype") +
  labs(x = "fastANI (%)", y = "Phenotype score",
       subtitle = sprintf(
         "Minimum observed fastANI. Red: %s (%.2f%%). Green: any kin (%.2f%%).",
         top_label, strong_kin_min, any_kin_min)) +
  theme_minimal(base_size = plot_base) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "threshold_plot.png"), p_thresh,
       width  = cfg$plots$scatter$width,
       height = cfg$plots$scatter$height, dpi = plot_dpi)

# --- 6c. Within- vs. between-group phenotype composition ----------------------
group_df <- data |>
  mutate(
    comparison = ifelse(Group_S1 == Group_S2, "Within group", "Between groups"),
    phenotype  = factor(
      phenotype_labels[as.character(Phenotype_interaction)],
      levels = phenotype_order
    )
  )

p_groups <- ggplot(group_df, aes(comparison, fill = phenotype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = setNames(unname(phenotype_colors),
                                      unname(phenotype_labels))) +
  labs(x = NULL, y = "Proportion of pairs", fill = "Phenotype") +
  theme_minimal(base_size = plot_base) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "within_vs_between_groups.png"), p_groups,
       width  = cfg$plots$groups$width,
       height = cfg$plots$groups$height, dpi = plot_dpi)

# --- 6d. fastANI distribution -------------------------------------------------
ani       <- data$fastANI[!is.na(data$fastANI)]
ani_stats <- tibble(mean   = mean(ani),
                    median = median(ani),
                    sd     = sd(ani),
                    min    = min(ani),
                    max    = max(ani),
                    n      = length(ani))

print(kable(ani_stats, digits = 2,
            caption = "fastANI distribution across all pairs."))

p_hist <- ggplot(tibble(ANI = ani), aes(ANI)) +
  geom_histogram(bins = plot_bins, fill = "skyblue",
                 color = "black", alpha = 0.7) +
  geom_vline(xintercept = mean(ani), linetype = "dashed", color = "red") +
  labs(x = "fastANI (%)", y = "Frequency") +
  theme_minimal(base_size = plot_base)

p_dens <- ggplot(tibble(ANI = ani), aes(ANI)) +
  geom_density(fill = "lightcoral", color = "black", alpha = 0.7) +
  geom_vline(xintercept = mean(ani),   linetype = "dashed", color = "red") +
  geom_vline(xintercept = median(ani), linetype = "dotted", color = "blue") +
  labs(x = "fastANI (%)", y = "Density") +
  theme_minimal(base_size = plot_base)

ggsave(file.path(output_dir, "fastani_distribution.png"),
       arrangeGrob(p_hist, p_dens, nrow = 1),
       width  = cfg$plots$distribution$width,
       height = cfg$plots$distribution$height, dpi = plot_dpi)

# --- 7. Summary report --------------------------------------------------------
phenotype_counts <- data |>
  count(Phenotype_interaction) |>
  mutate(label = phenotype_labels[as.character(Phenotype_interaction)]) |>
  arrange(desc(Phenotype_interaction))

summary_lines <- c(
  sprintf("Configuration:      %s", config_file),
  sprintf("Input:              %s", input_file),
  "",
  sprintf("Total pairs:        %d", nrow(data)),
  sprintf("Unique strains:     %d", length(unique(c(data$S1, data$S2)))),
  "",
  "Phenotype counts:",
  sprintf("  %-18s %d", phenotype_counts$label, phenotype_counts$n),
  "",
  sprintf("fastANI:            mean %.2f%% (SD %.2f%%)",
          ani_stats$mean, ani_stats$sd),
  sprintf("%s minimum: %.2f%% fastANI", top_label, strong_kin_min),
  sprintf("Any-kin minimum:    %.2f%% fastANI", any_kin_min),
  "  (minimum observed values; descriptive, not fitted. See S15B.)",
  "",
  "Nonconforming pairs:",
  sprintf("  %-48s %d", lab_strong,   nrow(strong_nonconformers)),
  sprintf("  %-48s %d", lab_moderate, nrow(moderate_nonconformers)),
  sprintf("  %-48s %d", lab_border,   nrow(borderline_pairs)),
  "",
  sprintf("Pearson r = %.4f (95%% CI %.4f to %.4f), p = %.2e",
          cor_res$estimate, cor_res$conf.int[1],
          cor_res$conf.int[2], cor_res$p.value)
)

writeLines(summary_lines, file.path(output_dir, "summary.txt"))
cat(summary_lines, sep = "\n")
cat("\n")
