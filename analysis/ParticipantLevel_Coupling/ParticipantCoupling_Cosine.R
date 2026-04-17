library(dplyr)
library(tidyr)

# Argument 1: cosine trial CSV. Argument 2: ventral trial CSV. Argument 3: output directory.
default_cosine <- "data/processed/gaze_target_cosine_occlusion.csv"
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_output_dir <- "outputs/Trial_decoupling"

args <- commandArgs(trailingOnly = TRUE)
cosine_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_cosine
ventral_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_ventral
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_csv <- file.path(output_dir, "within_phase_participant_coupling_cosine_summary.csv")
output_txt <- file.path(output_dir, "within_phase_participant_coupling_cosine.txt")

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

normalize_phase <- function(x) {
  case_when(
    x == "phase_i" ~ "implicit",
    x == "phase_e" ~ "explicit",
    x == "phase_p" ~ "perception",
    TRUE ~ as.character(x)
  )
}

add_perception_z <- function(df, dv_col) {
  df <- df %>%
    filter(phase %in% c("perception", "explicit", "implicit")) %>%
    filter(!is.na(.data[[dv_col]]), !is.na(participant), !is.na(stimulus))
  base <- df %>%
    filter(phase == "perception") %>%
    group_by(participant, stimulus) %>%
    summarize(
      base_mean = mean(.data[[dv_col]], na.rm = TRUE),
      base_sd = sd(.data[[dv_col]], na.rm = TRUE),
      .groups = "drop"
    )
  df %>%
    left_join(base, by = c("participant", "stimulus")) %>%
    mutate(dv_z = (.data[[dv_col]] - base_mean) / base_sd) %>%
    filter(is.finite(dv_z))
}

cosine <- read.csv(cosine_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = normalize_phase(tolower(phase)),
    stimulus = tolower(stimulus),
    cos_median = as.numeric(cos_median)
  ) %>%
  filter(is.finite(cos_median), !is.na(participant), !is.na(phase), !is.na(stimulus))

cosine_sum <- cosine %>%
  group_by(participant, phase) %>%
  summarize(dorsal_mean = mean(cos_median, na.rm = TRUE), .groups = "drop")

ventral <- read.csv(ventral_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    dv_raw = ifelse(stimulus == "white", -pupil_response, pupil_response)
  )

ventral_z <- add_perception_z(ventral, "dv_raw")

ventral_sum <- ventral_z %>%
  group_by(participant, phase) %>%
  summarize(ventral_mean = mean(dv_z, na.rm = TRUE), .groups = "drop")

summary_df <- full_join(cosine_sum, ventral_sum, by = c("participant", "phase")) %>%
  filter(phase %in% c("explicit", "implicit"))

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

phase_cor <- function(data, phase_value) {
  d <- data %>% filter(phase == phase_value)
  cor_stats(d$dorsal_mean, d$ventral_mean)
}

cor_exp <- phase_cor(summary_df, "explicit")
cor_imp <- phase_cor(summary_df, "implicit")

compare_bootstrap <- function(wide_df, n_boot = 5000, seed = 1) {
  d <- wide_df %>%
    filter(
      is.finite(dorsal_mean_explicit),
      is.finite(ventral_mean_explicit),
      is.finite(dorsal_mean_implicit),
      is.finite(ventral_mean_implicit)
    )
  n <- nrow(d)
  if (n < 3) return(list(n = n, r_exp = NA_real_, r_imp = NA_real_, delta = NA_real_, ci = c(NA, NA), p = NA_real_))

  x_exp <- d$dorsal_mean_explicit
  y_exp <- d$ventral_mean_explicit
  x_imp <- d$dorsal_mean_implicit
  y_imp <- d$ventral_mean_implicit

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

boot <- compare_bootstrap(wide)

fmt <- function(x, digits = 3) {
  if (!is.finite(x)) return("NA")
  formatC(x, digits = digits, format = "f")
}

summary_lines <- c(
  "=== Within-phase participant-level coupling (cosine dorsal) ===",
  "",
  "Mean-based (dorsal=cos_median, ventral=perception-z):",
  sprintf("  explicit: r = %s, p = %s, n = %d", fmt(cor_exp$r), fmt(cor_exp$p), cor_exp$n),
  sprintf("  implicit: r = %s, p = %s, n = %d", fmt(cor_imp$r), fmt(cor_imp$p), cor_imp$n),
  sprintf("  compare r_exp vs r_imp (bootstrap, n=%d): delta = %s, 95%% CI [%s, %s], p = %s",
          boot$n, fmt(boot$delta), fmt(boot$ci[1]), fmt(boot$ci[2]), fmt(boot$p))
)

writeLines(summary_lines, con = output_txt)
cat("Saved:", output_csv, "\n")
cat("Saved:", output_txt, "\n")
