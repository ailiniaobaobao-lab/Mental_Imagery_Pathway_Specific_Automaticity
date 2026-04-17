library(dplyr)
library(readxl)
library(ggplot2)

# Argument 1: gaze cosine CSV. Argument 2: participants xlsx. Argument 3: output directory.
default_csv <- "data/processed/gaze_target_cosine_occlusion.csv"
default_info <- "data/metadata/Participants_info.xlsx"
default_output_dir <- "outputs/Gaze_cos"

args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_csv
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Use participant-centered values instead of perception-z
use_participant_centered <- TRUE
metric_col <- "cos_mean"

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

info <- readxl::read_excel(info_path)
col_h <- names(info)[8]
group_map <- info %>%
  transmute(
    participant = pad_id(.data[["subject#"]]),
    group = as.character(.data[[col_h]])
  ) %>%
  filter(group %in% c("B", "N", "E"))

gaze_raw <- read.csv(csv_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = tolower(phase),
    stimulus = tolower(stimulus)
  ) %>%
  filter(phase %in% c("perception", "implicit", "explicit")) %>%
  inner_join(group_map, by = "participant")

metric_df <- gaze_raw %>%
  filter(!is.na(.data[[metric_col]]), n_samples > 0) %>%
  mutate(metric_value = .data[[metric_col]])

plot_df <- metric_df %>%
  group_by(participant) %>%
  mutate(
    participant_mean = mean(metric_value, na.rm = TRUE),
    dv = if (use_participant_centered) metric_value - participant_mean else metric_value
  ) %>%
  ungroup() %>%
  filter(is.finite(dv)) %>%
  group_by(group, phase) %>%
  summarize(
    mean = mean(dv, na.rm = TRUE),
    se = sd(dv, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    phase = factor(phase, levels = c("perception", "explicit", "implicit")),
    group = factor(group, levels = c("B", "E", "N"))
  )

p <- ggplot(plot_df, aes(x = phase, y = mean, color = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.08, linewidth = 0.6) +
  labs(
    title = "Gaze cosine mean by phase (participant-centered)",
    x = "Phase",
    y = "Cosine (centered)",
    color = "Group"
  ) +
  theme_minimal(base_size = 12)

out_path <- file.path(output_dir, "gaze_cosine_group_phase_mean.png")
ggsave(out_path, p, width = 6.5, height = 4.5, dpi = 300)
cat("Saved:", out_path, "\n")
