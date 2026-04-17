library(dplyr)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}
library(readxl)

# Argument 1: gaze cosine CSV. Argument 2: participants xlsx. Argument 3: output directory.
default_csv <- "data/processed/gaze_target_cosine_occlusion.csv"
default_info <- "data/metadata/Participants_info.xlsx"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_csv
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Use participant-centered values instead of perception-z
use_participant_centered <- TRUE

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

# Load participant group labels so cosine values can be modeled by subgroup.
info <- readxl::read_excel(info_path)
col_h <- names(info)[8]
group_map <- info %>%
  transmute(
    participant = pad_id(.data[["subject#"]]),
    group = as.character(.data[[col_h]])
  ) %>%
  filter(group %in% c("B", "N", "E"))

# Load trial-level cosine summaries and keep the phases used in the group model.
gaze_raw <- read.csv(csv_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = tolower(phase),
    stimulus = tolower(stimulus)
  ) %>%
  filter(phase %in% c("perception", "implicit", "explicit")) %>%
  inner_join(group_map, by = "participant")

# Fit one mixed model per cosine metric and export the model summary.
run_model <- function(metric_col, metric_label, output_file) {
  metric_df <- gaze_raw %>%
    filter(!is.na(.data[[metric_col]]), n_samples > 0) %>%
    mutate(metric_value = .data[[metric_col]])

  model_df <- metric_df %>%
    group_by(participant) %>%
    mutate(
      participant_mean = mean(metric_value, na.rm = TRUE),
      dv = if (use_participant_centered) metric_value - participant_mean else metric_value
    ) %>%
    ungroup() %>%
    filter(is.finite(dv)) %>%
    mutate(
      participant = factor(participant),
      phase = factor(phase, levels = c("explicit", "implicit", "perception")),
      group = factor(group, levels = c("B", "N", "E")),
      stimulus = factor(stimulus, levels = c("white", "black"))
    )

  model <- lmer(dv ~ phase * group + stimulus + (1 | participant), data = model_df, REML = FALSE)

  summary_text <- c(
    paste0(
      "=== Gaze cosine group mixed model: ", metric_label,
      if (use_participant_centered) " (participant-centered)" else " (raw)",
      " ~ phase * group + stimulus + (1|participant) ==="
    ),
    capture.output(summary(model))
  )

  writeLines(summary_text, con = output_file)
  cat("Saved:", output_file, "\n")
}

run_model(
  metric_col = "cos_mean",
  metric_label = "cos_mean",
  output_file = file.path(output_dir, "gaze_cosine_group_mixed_model_mean.txt")
)

run_model(
  metric_col = "cos_median",
  metric_label = "cos_median",
  output_file = file.path(output_dir, "gaze_cosine_group_mixed_model_median.txt")
)
