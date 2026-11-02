#!/usr/bin/env Rscript
# =============================================================================
# S15B_threshold_validation.R
#
# Uncertainty and predictive validation for the fastANI kin threshold.
#
# S15A reports the minimum observed fastANI per phenotype class and applies
# round-number cuts (99.4 / 99.5 / 99.7) to stratify nonconforming pairs.
# Those cuts are descriptive: a minimum is the most sample-dependent statistic
# available, and round numbers were chosen for readability, not fitted to data.
#
# This script quantifies how unstable that description is, and asks separately
# whether fastANI functions as a predictive boundary for kin phenotype:
#
#   1. Bootstrap CIs for the per-class minimum observed fastANI.
#      Resampling is done at the strain level, not the pair level, because
#      pairs sharing a strain are not independent observations.
#   2. Sensitivity sweep across 99.3-99.8% with sensitivity, specificity,
#      PPV, NPV and counts at every cut.
#   3. ROC (AUC + DeLong CI) and precision-recall analysis. PR is reported
#      alongside ROC because the kin/non-kin classes are imbalanced, which
#      inflates apparent ROC performance.
#   4. Leave-group-out validation: threshold selected on held-in groups
#      (Youden), evaluated on the held-out group.
#
# Inputs : data/strain_interactions.xlsx   (same file as S15A)
# Outputs: output/threshold_bootstrap_ci.csv
#          output/threshold_sensitivity_sweep.csv
#          output/threshold_roc_pr.csv
#          output/threshold_leave_group_out.csv
#          output/threshold_validation.png
#          output/threshold_validation_summary.txt
#
# Author: Eva Stare
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(knitr)
  library(gridExtra)
  library(pROC)
  library(tidyr)
  library(yaml)
})

# =============================================================================
# CONFIGURATION LAYER
# Paths, thresholds and settings are read from the shared configuration file
# used by S15A, so that the two scripts cannot drift apart.
#
# USAGE
#   Rscript S15B_threshold_validation.R [options]
#     --config=FILE     configuration file    (default: config_S15.yml)
#     --input=FILE      input .xlsx           (overrides config)
#     --output=DIR      output directory      (overrides config)
#     --kin-cutoff=X    phenotype >= X is kin (overrides config)
#     --n-boot=N        bootstrap replicates  (overrides config)
#     --help            show this message
# =============================================================================

usage <- function() {
  cat("
Usage: Rscript S15B_threshold_validation.R [options]

  --config=FILE     configuration file    (default: config_S15.yml)
  --input=FILE      input .xlsx           (overrides config)
  --output=DIR      output directory      (overrides config)
  --kin-cutoff=X    phenotype >= X is kin (overrides config)
  --n-boot=N        bootstrap replicates  (overrides config)
  --help            show this message

")
  quit(status = 0)
}

parse_args <- function(args) {
  if ("--help" %in% args || "-h" %in% args) usage()
  out <- list()
  for (a in args) {
    if (!grepl("^--[a-z-]+=", a))
      stop("Unrecognised argument: ", a,
           "\n  Run with --help for usage.", call. = FALSE)
    out[[sub("^--([a-z-]+)=.*$", "\\1", a)]] <- sub("^--[a-z-]+=", "", a)
  }
  out
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

config_file <- if (!is.null(opts$config)) opts$config else "config_S15.yml"
if (!file.exists(config_file))
  stop("Configuration file not found: ", config_file,
       "\n  Expected config_S15.yml in the working directory, or pass",
       " --config=FILE.", call. = FALSE)
cfg <- yaml::read_yaml(config_file)

input_file <- if (!is.null(opts$input))  opts$input  else cfg$paths$input_file
output_dir <- if (!is.null(opts$output)) opts$output else cfg$paths$output_dir

KIN_CUTOFF <- if (!is.null(opts[["kin-cutoff"]]))
                as.numeric(opts[["kin-cutoff"]]) else cfg$phenotype$kin_cutoff

N_BOOT     <- if (!is.null(opts[["n-boot"]]))
                as.integer(opts[["n-boot"]]) else cfg$validation$n_boot

SWEEP_RANGE <- seq(cfg$validation$sweep_from,
                   cfg$validation$sweep_to,
                   by = cfg$validation$sweep_by)

# Cuts used descriptively in S15A, carried through for comparison.
REFERENCE_CUTS <- as.numeric(cfg$validation$reference_cuts)

# Name of the pooled kin class, derived from KIN_CUTOFF so that plotting and
# summary code never references a hard-coded class label.
ANYKIN <- sprintf("Any kin (>= %.2g)", KIN_CUTOFF)

# Phenotype scale, expanded from the shared config.
phenotype_values <- vapply(cfg$phenotype$levels,
                           function(x) as.numeric(x$value), numeric(1))
phenotype_labels <- setNames(
  vapply(cfg$phenotype$levels, function(x) x$label, character(1)),
  as.character(phenotype_values))

set.seed(cfg$validation$seed)   # reproducible bootstrap

if (!file.exists(input_file))
  stop("Input file not found: ", input_file,
       "\n  Set paths:input_file in ", config_file,
       " or pass --input=FILE.", call. = FALSE)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Configuration:      ", config_file, "\n", sep = "")
cat("Input:              ", input_file,  "\n", sep = "")
cat("Output directory:   ", output_dir,  "\n", sep = "")


# --- 2. Load data -------------------------------------------------------------
data <- read_excel(input_file) |>
  mutate(
    Phenotype_interaction = as.numeric(Phenotype_interaction),
    fastANI               = as.numeric(fastANI)
  ) |>
  filter(!is.na(fastANI), !is.na(Phenotype_interaction)) |>
  mutate(is_kin = as.integer(Phenotype_interaction >= KIN_CUTOFF))

strains <- unique(c(data$S1, data$S2))

cat(sprintf("Pairs (complete cases): %d\n", nrow(data)))
cat(sprintf("Unique strains:         %d\n", length(strains)))
cat(sprintf("Kin pairs (>= %.2f):     %d (%.1f%%)\n",
            KIN_CUTOFF, sum(data$is_kin),
            100 * mean(data$is_kin)))
cat(sprintf("Non-kin pairs:          %d (%.1f%%)\n\n",
            sum(!data$is_kin), 100 * mean(!data$is_kin)))

if (length(unique(data$is_kin)) < 2)
  stop("Both kin and non-kin pairs are required for ROC/PR analysis.")

# --- 3. Bootstrap CIs for minimum observed fastANI ----------------------------
# Strain-level (cluster) bootstrap: resample strains with replacement, then
# keep pairs whose BOTH members are in the resampled set. This respects the
# dependence between pairs sharing a strain. Pair-level resampling would give
# CIs that are too narrow.

class_masks <- c(
  setNames(
    lapply(phenotype_values, function(v) {
      force(v); function(d) d$Phenotype_interaction == v
    }),
    unname(phenotype_labels[as.character(phenotype_values)])
  ),
  setNames(
    list(function(d) d$Phenotype_interaction >= KIN_CUTOFF),
    sprintf("Any kin (>= %.2g)", KIN_CUTOFF)
  )
)

boot_strain_sample <- function(d, all_strains) {
  keep <- sample(all_strains, length(all_strains), replace = TRUE)
  d[d$S1 %in% keep & d$S2 %in% keep, , drop = FALSE]
}

boot_mins <- matrix(NA_real_, nrow = N_BOOT, ncol = length(class_masks),
                    dimnames = list(NULL, names(class_masks)))

for (b in seq_len(N_BOOT)) {
  db <- boot_strain_sample(data, strains)
  if (!nrow(db)) next
  for (cl in names(class_masks)) {
    v <- db$fastANI[class_masks[[cl]](db)]
    if (length(v)) boot_mins[b, cl] <- min(v)
  }
}

bootstrap_ci <- bind_rows(lapply(names(class_masks), function(cl) {
  observed <- {
    v <- data$fastANI[class_masks[[cl]](data)]
    if (length(v)) min(v) else NA_real_
  }
  bv <- boot_mins[, cl]
  bv <- bv[!is.na(bv)]
  tibble(
    class          = cl,
    n_pairs        = sum(class_masks[[cl]](data)),
    observed_min   = observed,
    boot_median    = if (length(bv)) median(bv) else NA_real_,
    ci_lo_2.5      = if (length(bv)) unname(quantile(bv, 0.025)) else NA_real_,
    ci_hi_97.5     = if (length(bv)) unname(quantile(bv, 0.975)) else NA_real_,
    ci_width       = if (length(bv))
                       unname(diff(quantile(bv, c(0.025, 0.975)))) else NA_real_,
    n_boot_valid   = length(bv)
  )
}))

write.csv(bootstrap_ci,
          file.path(output_dir, "threshold_bootstrap_ci.csv"),
          row.names = FALSE)

print(kable(bootstrap_ci, digits = 3,
            caption = paste0("Bootstrap CIs for minimum observed fastANI ",
                             "(strain-level resampling, ", N_BOOT, " reps).")))

# --- 4. Sensitivity sweep across candidate thresholds -------------------------
# At each cut t, the rule "fastANI >= t predicts kin" is scored against the
# observed phenotype. PPV and NPV depend on the kin prevalence in this panel
# and will not transfer to a collection with different composition.

confusion_at <- function(d, t) {
  pred <- d$fastANI >= t
  tp <- sum( pred &  d$is_kin == 1)
  fp <- sum( pred &  d$is_kin == 0)
  fn <- sum(!pred &  d$is_kin == 1)
  tn <- sum(!pred &  d$is_kin == 0)
  tibble(
    threshold   = t,
    TP = tp, FP = fp, FN = fn, TN = tn,
    sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
    specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
    PPV         = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
    NPV         = if ((tn + fn) > 0) tn / (tn + fn) else NA_real_,
    accuracy    = (tp + tn) / nrow(d),
    youden      = (if ((tp + fn) > 0) tp / (tp + fn) else NA_real_) +
                  (if ((tn + fp) > 0) tn / (tn + fp) else NA_real_) - 1
  )
}

sweep <- bind_rows(lapply(SWEEP_RANGE, function(t) confusion_at(data, t))) |>
  mutate(used_in_S15 = threshold %in% REFERENCE_CUTS)

write.csv(sweep,
          file.path(output_dir, "threshold_sensitivity_sweep.csv"),
          row.names = FALSE)

print(kable(sweep, digits = 3,
            caption = "Sensitivity analysis across candidate fastANI cuts."))

best_youden_sweep <- sweep$threshold[which.max(sweep$youden)]

# --- 5. ROC and precision-recall ----------------------------------------------
roc_obj <- roc(response  = data$is_kin,
               predictor = data$fastANI,
               levels    = c(0, 1),
               direction = "<",
               quiet     = TRUE)

auc_ci <- ci.auc(roc_obj, method = "delong")

# Youden-optimal cut from the full ROC (not restricted to the sweep grid).
roc_best <- coords(roc_obj, "best", best.method = "youden",
                   ret = c("threshold", "sensitivity", "specificity"),
                   transpose = FALSE)
roc_best_threshold <- roc_best$threshold[1]

# Precision-recall computed directly from ROC coordinates.
prev <- mean(data$is_kin)
co <- coords(roc_obj, "all",
             ret = c("threshold", "sensitivity", "specificity"),
             transpose = FALSE) |>
  as_tibble() |>
  filter(is.finite(threshold)) |>
  arrange(threshold) |>
  mutate(
    recall    = sensitivity,
    precision = ifelse(
      (sensitivity * prev + (1 - specificity) * (1 - prev)) > 0,
      (sensitivity * prev) /
        (sensitivity * prev + (1 - specificity) * (1 - prev)),
      NA_real_
    )
  )

# Average precision (step-wise integral over recall).
pr_tbl <- co |> filter(!is.na(precision)) |> arrange(desc(recall))
avg_precision <- if (nrow(pr_tbl) > 1) {
  sum(diff(rev(pr_tbl$recall)) * rev(pr_tbl$precision)[-1])
} else NA_real_

roc_pr_summary <- tibble(
  metric = c("AUC (ROC)", "AUC CI lower (DeLong 95%)",
             "AUC CI upper (DeLong 95%)", "Average precision (PR)",
             "Baseline precision (kin prevalence)",
             "Youden-optimal fastANI cut"),
  value  = c(as.numeric(auc(roc_obj)), auc_ci[1], auc_ci[3],
             avg_precision, prev, roc_best_threshold)
)

write.csv(roc_pr_summary,
          file.path(output_dir, "threshold_roc_pr.csv"),
          row.names = FALSE)

print(kable(roc_pr_summary, digits = 4,
            caption = "ROC and precision-recall summary."))

# --- 6. Leave-group-out validation --------------------------------------------
# Held-out unit = candidate KD group (Group_S1/Group_S2 in the input table)

data <- data |> mutate(unit_S1 = Group_S1, unit_S2 = Group_S2)
unit_label <- "candidate KD group"

units <- sort(unique(na.omit(c(data$unit_S1, data$unit_S2))))

lgo <- bind_rows(lapply(units, function(u) {
  test  <- data |> filter(unit_S1 == u | unit_S2 == u)
  train <- data |> filter(unit_S1 != u,  unit_S2 != u)

  if (!nrow(test) || !nrow(train) ||
      length(unique(train$is_kin)) < 2)
    return(tibble(held_out = as.character(u),
                  n_train = nrow(train), n_test = nrow(test),
                  fitted_threshold = NA_real_,
                  test_sensitivity = NA_real_, test_specificity = NA_real_,
                  test_PPV = NA_real_, test_NPV = NA_real_,
                  test_accuracy = NA_real_,
                  note = "insufficient data / single class in training fold"))

  r_train <- roc(train$is_kin, train$fastANI,
                 levels = c(0, 1), direction = "<", quiet = TRUE)
  t_fit <- coords(r_train, "best", best.method = "youden",
                  ret = "threshold", transpose = FALSE)$threshold[1]

  cm <- confusion_at(test, t_fit)
  tibble(
    held_out         = as.character(u),
    n_train          = nrow(train),
    n_test           = nrow(test),
    fitted_threshold = t_fit,
    test_sensitivity = cm$sensitivity,
    test_specificity = cm$specificity,
    test_PPV         = cm$PPV,
    test_NPV         = cm$NPV,
    test_accuracy    = cm$accuracy,
    note             = if (length(unique(test$is_kin)) < 2)
                         "test fold single-class; some metrics undefined" else ""
  )
}))

write.csv(lgo,
          file.path(output_dir, "threshold_leave_group_out.csv"),
          row.names = FALSE)

print(kable(lgo, digits = 3,
            caption = sprintf("Leave-group-out validation, held-out unit = %s.",
                              unit_label)))

thr_spread <- lgo$fitted_threshold[!is.na(lgo$fitted_threshold)]

# --- 7. Plots -----------------------------------------------------------------

# 7a. Bootstrap distribution of the "any kin" minimum.
boot_anykin <- boot_mins[, ANYKIN]
boot_anykin <- boot_anykin[!is.na(boot_anykin)]

p_boot <- ggplot(tibble(x = boot_anykin), aes(x)) +
  geom_histogram(bins = 40, fill = "#90A4AE", color = "black", alpha = 0.8) +
  geom_vline(xintercept = bootstrap_ci$observed_min[
               bootstrap_ci$class == ANYKIN],
             color = "darkgreen", linewidth = 1) +
  geom_vline(xintercept = bootstrap_ci$ci_lo_2.5[
               bootstrap_ci$class == ANYKIN],
             linetype = "dashed", color = "red") +
  geom_vline(xintercept = bootstrap_ci$ci_hi_97.5[
               bootstrap_ci$class == ANYKIN],
             linetype = "dashed", color = "red") +
  labs(title = "Bootstrap distribution: minimum fastANI among kin pairs",
       subtitle = "Green: observed minimum. Red dashed: 95% percentile interval.",
       x = "Minimum fastANI (%)", y = "Bootstrap replicates") +
  theme_minimal(base_size = 11)

# 7b. Sensitivity sweep.
sweep_long <- sweep |>
  select(threshold, sensitivity, specificity, PPV, NPV) |>
  tidyr::pivot_longer(-threshold, names_to = "metric", values_to = "value")

p_sweep <- ggplot(sweep_long, aes(threshold, value, color = metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  geom_vline(xintercept = REFERENCE_CUTS, linetype = "dotted",
             color = "grey40") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Performance across candidate fastANI thresholds",
       subtitle = "Grey dotted: cuts used descriptively in S15A.",
       x = "fastANI threshold (%)", y = "Value", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# 7c. ROC.
roc_df <- tibble(fpr = 1 - roc_obj$specificities,
                 tpr = roc_obj$sensitivities) |>
  arrange(fpr)

p_roc <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_line(linewidth = 0.9, color = "#1565C0") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60") +
  labs(title = sprintf("ROC (AUC = %.3f, 95%% CI %.3f-%.3f)",
                       as.numeric(auc(roc_obj)), auc_ci[1], auc_ci[3]),
       x = "False positive rate", y = "True positive rate") +
  theme_minimal(base_size = 11)

# 7d. Precision-recall.
p_pr <- ggplot(pr_tbl, aes(recall, precision)) +
  geom_line(linewidth = 0.9, color = "#C62828") +
  geom_hline(yintercept = prev, linetype = "dashed", color = "grey60") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = sprintf("Precision-recall (AP = %.3f)", avg_precision),
       subtitle = "Grey dashed: baseline precision = kin prevalence.",
       x = "Recall", y = "Precision") +
  theme_minimal(base_size = 11)

ggsave(file.path(output_dir, "threshold_validation.png"),
       arrangeGrob(p_boot, p_sweep, p_roc, p_pr, nrow = 2),
       width = 13, height = 9, dpi = 300)

# --- 8. Summary report --------------------------------------------------------
anykin_row <- bootstrap_ci |> filter(class == ANYKIN)

summary_lines <- c(
  "S15B: fastANI threshold uncertainty and predictive validation",
  strrep("=", 62),
  "",
  sprintf("Pairs analysed:       %d", nrow(data)),
  sprintf("Unique strains:       %d", length(strains)),
  sprintf("Kin prevalence:       %.1f%% (phenotype >= %.2f)",
          100 * prev, KIN_CUTOFF),
  "",
  "1. BOOTSTRAP CI FOR MINIMUM OBSERVED fastANI",
  sprintf("   (strain-level resampling, %d replicates)", N_BOOT),
  sprintf("   Any kin: observed min %.3f%%, 95%% CI %.3f-%.3f%% (width %.3f)",
          anykin_row$observed_min, anykin_row$ci_lo_2.5,
          anykin_row$ci_hi_97.5, anykin_row$ci_width),
  "",
  sprintf("2. SENSITIVITY ANALYSIS (%.2f-%.2f%%)",
          min(SWEEP_RANGE), max(SWEEP_RANGE)),
  sprintf("   Youden-optimal cut on the sweep grid: %.2f%%", best_youden_sweep),
  # Reference cuts come from the config; matched by nearest grid value rather
  # than floating-point equality, which fails if the sweep step changes.
  unlist(lapply(REFERENCE_CUTS, function(cut) {
    i <- which.min(abs(sweep$threshold - cut))
    if (abs(sweep$threshold[i] - cut) > 1e-8) return(NULL)
    sprintf("   At %.2f%%: sens %.3f, spec %.3f, PPV %.3f, NPV %.3f",
            sweep$threshold[i], sweep$sensitivity[i], sweep$specificity[i],
            sweep$PPV[i], sweep$NPV[i])
  })),
  "",
  "3. ROC AND PRECISION-RECALL",
  sprintf("   AUC = %.3f (95%% CI %.3f-%.3f, DeLong)",
          as.numeric(auc(roc_obj)), auc_ci[1], auc_ci[3]),
  sprintf("   Average precision = %.3f (baseline = %.3f)",
          avg_precision, prev),
  sprintf("   Youden-optimal cut (full ROC): %.3f%%", roc_best_threshold),
  "",
  sprintf("4. LEAVE-GROUP-OUT VALIDATION (held-out unit: %s)", unit_label),
  sprintf("   Folds evaluated: %d of %d units", length(thr_spread),
          length(units)),
  if (length(thr_spread))
    sprintf("   Fitted threshold across folds: %.3f-%.3f%% (median %.3f%%)",
            min(thr_spread), max(thr_spread), median(thr_spread))
  else "   No fold had sufficient data for threshold fitting.",
  sprintf("   Mean held-out accuracy: %.3f",
          mean(lgo$test_accuracy, na.rm = TRUE)),
  "   Candidate KD groups are used as the held-out unit because the interaction",
  "   table carries no lineage assignment."
)

writeLines(summary_lines,
           file.path(output_dir, "threshold_validation_summary.txt"))
cat("\n")
cat(summary_lines, sep = "\n")
cat("\n")