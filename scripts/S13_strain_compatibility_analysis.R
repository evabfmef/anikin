#!/usr/bin/env Rscript
# =============================================================================
# S13_strain_compatibility_analysis.R
#
# Kin compatibility: fastANI thresholds and within-group concordance.
#
# Analyses the concordance between pairwise swarming phenotypes and fastANI
# genomic similarity across the 78-strain panel, produces the fastANI
# threshold summary reported in the main text, generates the fastANI
# distribution plot (Figure 2C), and stratifies nonconforming pairs into
# strong / moderate / borderline classes.
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(knitr)
  library(gridExtra)
})

# --- 1. Parameters (edit as needed) -------------------------------------------
input_file <- "data/strain_interactions.xlsx"
output_dir <- "output"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

phenotype_colors <- c("1"    = "#4CAF50", "0.5"  = "#9CCC65",
                      "0.25" = "#FF9800", "0"    = "#F44336")
phenotype_labels <- c("1"    = "Strong kin", "0.5"  = "Kin",
                      "0.25" = "Non-kin",    "0"    = "Strong non-kin")

# --- 2. Load data -------------------------------------------------------------
data <- read_excel(input_file) |>
  mutate(
    Phenotype_interaction = as.numeric(Phenotype_interaction),
    fastANI               = as.numeric(fastANI)
  )

cat("Pairs loaded:",    nrow(data),                            "\n")
cat("Unique strains:",  length(unique(c(data$S1, data$S2))),   "\n")

# --- 3. Within-group concordance ----------------------------------------------
# For each KD group with >=2 tested members, summarise within-group swarming
# phenotypes and mean fastANI. Kin-consistent: all phenotypes >= 0.5.
# Weak-kin-consistent: all phenotypes >= 0.25 but some below 0.5.
within_group <- data |>
  filter(Group_S1 == Group_S2) |>
  group_by(group = Group_S1) |>
  summarise(
    n_pairs             = n(),
    min_phenotype       = min( Phenotype_interaction, na.rm = TRUE),
    mean_phenotype      = mean(Phenotype_interaction, na.rm = TRUE),
    mean_fastani        = mean(fastANI,               na.rm = TRUE),
    kin_consistent      = all(Phenotype_interaction >= 0.5,  na.rm = TRUE),
    weak_kin_consistent = all(Phenotype_interaction >= 0.25, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(mean_phenotype))

write.csv(within_group,
          file.path(output_dir, "within_group_summary.csv"),
          row.names = FALSE)

print(kable(within_group, digits = 2,
            caption = "Within-group swarming phenotype and fastANI summary."))

# --- 4. fastANI thresholds per phenotype class --------------------------------
# Minimum fastANI at which each phenotype class is observed defines the
# empirical threshold below which that phenotype is not seen in this dataset.
threshold_summary <- function(d, mask, label) {
  vals <- d$fastANI[mask & !is.na(d$fastANI)]
  if (!length(vals)) return(NULL)
  tibble(class = label, n = length(vals),
         min_fastani = min(vals), max_fastani = max(vals))
}

thr <- bind_rows(
  threshold_summary(data, data$Phenotype_interaction == 1,    "Strong kin"),
  threshold_summary(data, data$Phenotype_interaction == 0.5,  "Kin"),
  threshold_summary(data, data$Phenotype_interaction >= 0.5,  "Any kin (>= 0.5)"),
  threshold_summary(data, data$Phenotype_interaction == 0.25, "Non-kin"),
  threshold_summary(data, data$Phenotype_interaction == 0,    "Strong non-kin")
)

print(kable(thr, digits = 2,
            caption = "fastANI ranges by phenotype class."))

# --- 5. Pearson correlation ---------------------------------------------------
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

# --- 6. Nonconforming pairs stratification ------------------------------------
# Strong nonconformers    : fastANI >= 99.7% AND phenotype = 0 (lysis)
# Moderate nonconformers  : fastANI >= 99.5% AND phenotype <= 0.25
# Borderline pairs        : 99.4% <= fastANI < 99.5%
strong_nonconformers <- data |>
  filter(fastANI >= 99.7, Phenotype_interaction == 0) |>
  arrange(desc(fastANI))

moderate_nonconformers <- data |>
  filter(fastANI >= 99.5, Phenotype_interaction <= 0.25) |>
  arrange(desc(fastANI), Phenotype_interaction)

borderline_pairs <- data |>
  filter(fastANI >= 99.4, fastANI < 99.5) |>
  arrange(desc(fastANI))

nonconf_counts <- tibble(
  Category = c("Strong nonconformers (>=99.7% + lysis)",
               "Moderate nonconformers (>=99.5% + any non-kin)",
               "Borderline pairs (99.4-99.5% ANI)"),
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
      kable(digits = 3,
            caption = "Strong nonconforming pairs (fastANI >= 99.7% + lysis).")
  )
}

# --- 7. Plots -----------------------------------------------------------------

# --- 7a. fastANI vs. phenotype ------------------------------------------------
p_scatter <- ggplot(data, aes(fastANI, Phenotype_interaction)) +
  geom_jitter(aes(color = factor(Phenotype_interaction)),
              size = 2.5, alpha = 0.7, height = 0.02, width = 0) +
  geom_smooth(method = "loess", se = TRUE, color = "black", alpha = 0.2) +
  scale_color_manual(values = phenotype_colors,
                     labels = phenotype_labels,
                     name   = "Phenotype") +
  labs(x = "fastANI (%)", y = "Phenotype score") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fastani_vs_phenotype.png"),
       p_scatter, width = 8, height = 5, dpi = 300)

# --- 7b. Threshold plot -------------------------------------------------------
strong_kin_min <- thr$min_fastani[thr$class == "Strong kin"]
any_kin_min    <- thr$min_fastani[thr$class == "Any kin (>= 0.5)"]

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
         "Red: strong-kin minimum (%.2f%%). Green: any-kin minimum (%.2f%%).",
         strong_kin_min, any_kin_min)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "threshold_plot.png"),
       p_thresh, width = 8, height = 5, dpi = 300)

# --- 7c. Within- vs. between-group phenotype composition ----------------------
group_df <- data |>
  mutate(
    comparison = ifelse(Group_S1 == Group_S2, "Within group", "Between groups"),
    phenotype  = factor(
      phenotype_labels[as.character(Phenotype_interaction)],
      levels = c("Strong kin", "Kin", "Non-kin", "Strong non-kin")
    )
  )

p_groups <- ggplot(group_df, aes(comparison, fill = phenotype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = setNames(phenotype_colors, phenotype_labels)) +
  labs(x = NULL, y = "Proportion of pairs", fill = "Phenotype") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "within_vs_between_groups.png"),
       p_groups, width = 7, height = 5, dpi = 300)

# --- 7d. fastANI distribution -------------------------------------------------
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
  geom_histogram(bins = 50, fill = "skyblue",
                 color = "black", alpha = 0.7) +
  geom_vline(xintercept = mean(ani), linetype = "dashed", color = "red") +
  labs(x = "fastANI (%)", y = "Frequency") +
  theme_minimal(base_size = 12)

p_dens <- ggplot(tibble(ANI = ani), aes(ANI)) +
  geom_density(fill = "lightcoral", color = "black", alpha = 0.7) +
  geom_vline(xintercept = mean(ani),   linetype = "dashed", color = "red") +
  geom_vline(xintercept = median(ani), linetype = "dotted", color = "blue") +
  labs(x = "fastANI (%)", y = "Density") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "fastani_distribution.png"),
       arrangeGrob(p_hist, p_dens, nrow = 1),
       width = 12, height = 4, dpi = 300)

# --- 8. Summary report --------------------------------------------------------
phenotype_counts <- data |>
  count(Phenotype_interaction) |>
  mutate(label = phenotype_labels[as.character(Phenotype_interaction)]) |>
  arrange(desc(Phenotype_interaction))

summary_lines <- c(
  sprintf("Total pairs:        %d", nrow(data)),
  sprintf("Unique strains:     %d", length(unique(c(data$S1, data$S2)))),
  "",
  "Phenotype counts:",
  sprintf("  %-18s %d", phenotype_counts$label, phenotype_counts$n),
  "",
  sprintf("fastANI:            mean %.2f%% (SD %.2f%%)",
          ani_stats$mean, ani_stats$sd),
  sprintf("Strong-kin minimum: %.2f%% fastANI", strong_kin_min),
  sprintf("Any-kin minimum:    %.2f%% fastANI", any_kin_min),
  "",
  "Nonconforming pairs:",
  sprintf("  Strong (>=99.7%% + lysis):          %d",
          nrow(strong_nonconformers)),
  sprintf("  Moderate (>=99.5%% + any non-kin):  %d",
          nrow(moderate_nonconformers)),
  sprintf("  Borderline (99.4-99.5%% ANI):       %d",
          nrow(borderline_pairs)),
  "",
  sprintf("Pearson r = %.4f (95%% CI %.4f to %.4f), p = %.2e",
          cor_res$estimate, cor_res$conf.int[1],
          cor_res$conf.int[2], cor_res$p.value)
)

writeLines(summary_lines, file.path(output_dir, "summary.txt"))
cat(summary_lines, sep = "\n")
