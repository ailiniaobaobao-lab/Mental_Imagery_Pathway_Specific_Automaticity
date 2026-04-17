library(dplyr)
library(tidyr)
library(lme4)

# ==== Input/Output ====
# Argument 1: ventral trial CSV. Argument 2: cosine trial CSV. Argument 3: optional output directory.
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_dorsal <- "data/processed/gaze_target_cosine_occlusion.csv"
default_output_dir <- "outputs/Trial_decoupling"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
cosine_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_dorsal
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_trials <- file.path(output_dir, "automaticity_trial_coupling_cosine_data.csv")
output_txt <- file.path(output_dir, "automaticity_trial_coupling_cosine.txt")

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

to_string_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return(sprintf("%.2e", p))
  sprintf("%.4f", p)
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

# ==== Ventral: white sign-flip + perception standardization ====
ventral <- read.csv(ventral_path, stringsAsFactors = FALSE)
ventral_clean <- ventral %>%
  mutate(
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    participant = pad_id(participant),
    trial_id = sub(".*-(\\d+)\\..*$", "\\1", as.character(video_file_mapped)),
    dv_raw = ifelse(stimulus == "white", -pupil_response, pupil_response)
  )

ventral_z <- add_perception_z(ventral_clean, "dv_raw") %>%
  filter(!is.na(trial_id)) %>%
  filter(phase %in% c("explicit", "implicit")) %>%
  group_by(participant, phase, stimulus, trial_id) %>%
  summarize(ventral_z = mean(dv_z, na.rm = TRUE), .groups = "drop")

# ==== Dorsal (cosine): cos_median (raw) ====
cosine <- read.csv(cosine_path, stringsAsFactors = FALSE)
cosine_clean <- cosine %>%
  mutate(
    phase = normalize_phase(tolower(phase)),
    stimulus = tolower(stimulus),
    participant = pad_id(participant),
    trial_id = as.character(trial_id),
    cos_median = as.numeric(cos_median)
  )

cosine_z <- cosine_clean %>%
  filter(phase %in% c("explicit", "implicit")) %>%
  group_by(participant, phase, stimulus, trial_id) %>%
  summarize(cosine_z = mean(cos_median, na.rm = TRUE), .groups = "drop")

# ==== Trial-level Coupling ====
trial_df <- inner_join(
  ventral_z,
  cosine_z,
  by = c("participant", "phase", "stimulus", "trial_id")
) %>%
  mutate(
    phase = factor(phase, levels = c("explicit", "implicit")),
    stimulus = factor(stimulus, levels = c("white", "black"))
  ) %>%
  group_by(participant) %>%
  mutate(
    cosine_z_w = cosine_z - mean(cosine_z, na.rm = TRUE),
    ventral_z_w = ventral_z - mean(ventral_z, na.rm = TRUE)
  ) %>%
  ungroup()

write.csv(trial_df, output_trials, row.names = FALSE)

summary_lines <- c()

if (requireNamespace("rmcorr", quietly = TRUE)) {
  library(rmcorr)
  rmc <- try(rmcorr(participant, cosine_z, ventral_z, trial_df), silent = TRUE)
  if (!inherits(rmc, "try-error")) {
    summary_lines <- c(
      "=== rmcorr (within-subject coupling) ===",
      paste0(
        "r = ", format(rmc$r, digits = 4),
        ", 95% CI [", format(rmc$CI[1], digits = 4),
        ", ", format(rmc$CI[2], digits = 4),
        "], p = ", to_string_p(rmc$p)
      )
    )
  } else {
    summary_lines <- c("=== rmcorr failed (insufficient levels) ===")
  }
} else {
  summary_lines <- c("=== rmcorr not available (package 'rmcorr' not installed) ===")
}

# Overall centered Pearson r (auxiliary)
overall_r <- stats::cor(trial_df$cosine_z_w, trial_df$ventral_z_w, use = "complete.obs")
overall_n <- sum(is.finite(trial_df$cosine_z_w) & is.finite(trial_df$ventral_z_w))
summary_lines <- c(
  summary_lines,
  "",
  "=== Auxiliary centered Pearson r (overall) ===",
  paste0("r = ", format(overall_r, digits = 4), ", n = ", overall_n)
)

model <- lmer(ventral_z_w ~ cosine_z_w + (1 | participant), data = trial_df, REML = FALSE)
ci_method <- "Wald"
ci <- try(confint(model, method = ci_method), silent = TRUE)
coef_df <- as.data.frame(summary(model)$coefficients)
est <- if ("cosine_z_w" %in% rownames(coef_df)) coef_df["cosine_z_w", "Estimate"] else NA_real_
if (!inherits(ci, "try-error") && "cosine_z_w" %in% rownames(ci)) {
  ci_line <- paste0(
    "beta(cosine_z_w) = ", format(est, digits = 4),
    ", 95% CI (", ci_method, "): [",
    format(ci["cosine_z_w", 1], digits = 4),
    ", ",
    format(ci["cosine_z_w", 2], digits = 4),
    "]"
  )
} else {
  if ("cosine_z_w" %in% rownames(coef_df)) {
    se <- coef_df["cosine_z_w", "Std. Error"]
    ci_low <- est - 1.96 * se
    ci_high <- est + 1.96 * se
    ci_line <- paste0(
      "beta(cosine_z_w) = ", format(est, digits = 4),
      ", 95% CI (Wald): [",
      format(ci_low, digits = 4),
      ", ",
      format(ci_high, digits = 4),
      "]"
    )
  } else {
    ci_line <- "beta(cosine_z_w) 95% CI: NA"
  }
}
summary_lines <- c(
  summary_lines,
  "",
  "=== Mixed model: ventral_z_w ~ cosine_z_w + (1|participant) ===",
  capture.output(summary(model)),
  "",
  ci_line
)

# Overall model with phase + stimulus controls
model_ctrl <- lmer(ventral_z_w ~ cosine_z_w + phase + stimulus + (1 | participant), data = trial_df, REML = FALSE)
summary_lines <- c(
  summary_lines,
  "",
  "=== Mixed model (controls): ventral_z_w ~ cosine_z_w + phase + stimulus + (1|participant) ===",
  capture.output(summary(model_ctrl))
)

# Refit the coupling analyses within a single phase.
add_phase_specific <- function(phase_label) {
  sub_df <- trial_df %>%
    filter(phase == phase_label) %>%
    group_by(participant) %>%
    mutate(
      cosine_z_w = cosine_z - mean(cosine_z, na.rm = TRUE),
      ventral_z_w = ventral_z - mean(ventral_z, na.rm = TRUE)
    ) %>%
    ungroup()
  if (!nrow(sub_df)) return(NULL)

  lines <- c("", paste0("=== Phase-specific: ", phase_label, " ==="))

  if (requireNamespace("rmcorr", quietly = TRUE)) {
    rmc <- try(rmcorr(participant, cosine_z, ventral_z, sub_df), silent = TRUE)
    if (!inherits(rmc, "try-error")) {
      lines <- c(
        lines,
        "=== rmcorr (within-subject coupling) ===",
        paste0(
          "r = ", format(rmc$r, digits = 4),
          ", 95% CI [", format(rmc$CI[1], digits = 4),
          ", ", format(rmc$CI[2], digits = 4),
          "], p = ", to_string_p(rmc$p)
        )
      )
    } else {
      lines <- c(lines, "=== rmcorr failed (insufficient levels) ===")
    }
  } else {
    lines <- c(lines, "=== rmcorr not available (package 'rmcorr' not installed) ===")
  }

  phase_r <- stats::cor(sub_df$cosine_z_w, sub_df$ventral_z_w, use = "complete.obs")
  phase_n <- sum(is.finite(sub_df$cosine_z_w) & is.finite(sub_df$ventral_z_w))
  lines <- c(
    lines,
    "",
    "=== Auxiliary centered Pearson r ===",
    paste0("r = ", format(phase_r, digits = 4), ", n = ", phase_n)
  )

  phase_model <- lmer(ventral_z_w ~ cosine_z_w + (1 | participant), data = sub_df, REML = FALSE)
  ci <- try(confint(phase_model, method = ci_method), silent = TRUE)
  coef_df <- as.data.frame(summary(phase_model)$coefficients)
  est <- if ("cosine_z_w" %in% rownames(coef_df)) coef_df["cosine_z_w", "Estimate"] else NA_real_
  if (!inherits(ci, "try-error") && "cosine_z_w" %in% rownames(ci)) {
    ci_line <- paste0(
      "beta(cosine_z_w) = ", format(est, digits = 4),
      ", 95% CI (", ci_method, "): [",
      format(ci["cosine_z_w", 1], digits = 4),
      ", ",
      format(ci["cosine_z_w", 2], digits = 4),
      "]"
    )
  } else {
    if ("cosine_z_w" %in% rownames(coef_df)) {
      se <- coef_df["cosine_z_w", "Std. Error"]
      ci_low <- est - 1.96 * se
      ci_high <- est + 1.96 * se
      ci_line <- paste0(
        "beta(cosine_z_w) = ", format(est, digits = 4),
        ", 95% CI (Wald): [",
        format(ci_low, digits = 4),
        ", ",
        format(ci_high, digits = 4),
        "]"
      )
    } else {
      ci_line <- "beta(cosine_z_w) 95% CI: NA"
    }
  }

  lines <- c(
    lines,
    "",
    "=== Mixed model: ventral_z_w ~ cosine_z_w + (1|participant) ===",
    capture.output(summary(phase_model)),
    "",
    ci_line
  )
  lines
}

summary_lines <- c(
  summary_lines,
  add_phase_specific("explicit"),
  add_phase_specific("implicit")
)

writeLines(summary_lines, con = output_txt)
cat("Saved:", output_trials, "\n")
cat("Saved:", output_txt, "\n")
