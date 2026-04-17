library(dplyr)
library(tidyr)

default_input <- "outputs/Trial_decoupling/automaticity_trial_coupling_data.csv"
default_output_dir <- "outputs/Trial_decoupling"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_dir <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_csv <- file.path(output_dir, "within_phase_participant_coupling_summary.csv")
output_txt <- file.path(output_dir, "within_phase_participant_coupling.txt")

# Load trial-level coupling data and keep the phases compared in the summary.
df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = tolower(phase),
    participant = as.character(participant)
  ) %>%
  filter(phase %in% c("explicit", "implicit"))

# Collapse trial-level values to participant-level pathway means within each phase.
summary_df <- df %>%
  group_by(participant, phase) %>%
  summarize(
    dorsal_mean = mean(corridor_z, na.rm = TRUE),
    ventral_mean = mean(ventral_z, na.rm = TRUE),
    n_trials = sum(is.finite(corridor_z) & is.finite(ventral_z)),
    .groups = "drop"
  )

write.csv(summary_df, output_csv, row.names = FALSE)

cor_stats <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < 3) return(list(r = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  list(r = unname(ct$estimate), p = ct$p.value, n = n)
}

phase_cor <- function(data, dorsal_col, ventral_col, phase_value) {
  d <- data %>% filter(phase == phase_value)
  cor_stats(d[[dorsal_col]], d[[ventral_col]])
}

# Estimate the phase-specific correlations and compare them with bootstrap resampling.
cor_exp_mean <- phase_cor(summary_df, "dorsal_mean", "ventral_mean", "explicit")
cor_imp_mean <- phase_cor(summary_df, "dorsal_mean", "ventral_mean", "implicit")
compare_bootstrap <- function(wide_df, dorsal_prefix, ventral_prefix, n_boot = 5000, seed = 1) {
  d <- wide_df %>%
    filter(
      is.finite(.data[[paste0(dorsal_prefix, "_explicit")]]),
      is.finite(.data[[paste0(ventral_prefix, "_explicit")]]),
      is.finite(.data[[paste0(dorsal_prefix, "_implicit")]]),
      is.finite(.data[[paste0(ventral_prefix, "_implicit")]])
    )
  n <- nrow(d)
  if (n < 3) return(list(n = n, r_exp = NA_real_, r_imp = NA_real_, delta = NA_real_, ci = c(NA, NA), p = NA_real_))

  x_exp <- d[[paste0(dorsal_prefix, "_explicit")]]
  y_exp <- d[[paste0(ventral_prefix, "_explicit")]]
  x_imp <- d[[paste0(dorsal_prefix, "_implicit")]]
  y_imp <- d[[paste0(ventral_prefix, "_implicit")]]

  r_exp <- suppressWarnings(cor(x_exp, y_exp, use = "complete.obs"))
  r_imp <- suppressWarnings(cor(x_imp, y_imp, use = "complete.obs"))
  delta_obs <- r_exp - r_imp

  set.seed(seed)
  deltas <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    rx <- suppressWarnings(cor(x_exp[idx], y_exp[idx], use = "complete.obs"))
    ri <- suppressWarnings(cor(x_imp[idx], y_imp[idx], use = "complete.obs"))
    deltas[i] <- rx - ri
  }

  ci <- quantile(deltas, probs = c(0.025, 0.975), na.rm = TRUE)
  p_raw <- 2 * min(mean(deltas >= 0, na.rm = TRUE), mean(deltas <= 0, na.rm = TRUE))
  p <- min(1, p_raw)
  list(n = n, r_exp = r_exp, r_imp = r_imp, delta = delta_obs, ci = ci, p = p)
}

wide <- summary_df %>%
  pivot_wider(
    names_from = phase,
    values_from = c(dorsal_mean, ventral_mean),
    names_sep = "_"
  )

boot_mean <- compare_bootstrap(wide, "dorsal_mean", "ventral_mean")

fmt <- function(x, digits = 3) {
  if (!is.finite(x)) return("NA")
  formatC(x, digits = digits, format = "f")
}

summary_lines <- c(
  "=== Within-phase participant-level coupling ===",
  "",
  "Mean-based (dorsal_mean vs ventral_mean):",
  sprintf("  explicit: r = %s, p = %s, n = %d", fmt(cor_exp_mean$r), fmt(cor_exp_mean$p), cor_exp_mean$n),
  sprintf("  implicit: r = %s, p = %s, n = %d", fmt(cor_imp_mean$r), fmt(cor_imp_mean$p), cor_imp_mean$n),
  sprintf("  compare r_exp vs r_imp (bootstrap, n=%d): delta = %s, 95%% CI [%s, %s], p = %s",
          boot_mean$n, fmt(boot_mean$delta), fmt(boot_mean$ci[1]), fmt(boot_mean$ci[2]), fmt(boot_mean$p))
)

writeLines(summary_lines, con = output_txt)
cat("Saved:", output_csv, "\n")
cat("Saved:", output_txt, "\n")
