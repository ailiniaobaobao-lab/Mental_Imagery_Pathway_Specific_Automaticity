library(dplyr)
library(ggplot2)
library(tidyr)

# ==== Input/Output ====
default_input <- "data/processed/RDS_Data"
default_output_dir <- "outputs/PupilSize_Raw/PupilSize_mean_ci"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_dir <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

bin_ms <- 40
x_limit <- 30000
x_breaks <- seq(0, x_limit, by = 5000)

phase_config <- list(
  perception = list(phase_tag = "phase_p", video_col = "video_file_p", trial_col = "TRIAL_INDEX", suffix_keep = "0"),
  implicit   = list(phase_tag = "phase_i", video_col = "video_file_i", trial_col = "TRIAL_INDEX", suffix_keep = "1"),
  explicit   = list(phase_tag = "phase_e", video_col = "video_file_e", trial_col = "TRIAL_INDEX", suffix_keep = "1")
)

parse_video_info <- function(video_name) {
  base <- tools::file_path_sans_ext(basename(video_name))
  stimulus <- if (grepl("white", base, ignore.case = TRUE)) {
    "white"
  } else if (grepl("black", base, ignore.case = TRUE)) {
    "black"
  } else {
    NA_character_
  }
  suffix <- sub(".*-(\\d+)$", "\\1", base)
  if (identical(suffix, base)) suffix <- NA_character_
  list(stimulus = stimulus, suffix = suffix)
}

select_primary_events <- function(df, trial_col) {
  trial_sym <- rlang::sym(trial_col)
  df %>%
    group_by(!!trial_sym, Event) %>%
    summarize(
      samples = n(),
      trial_index = suppressWarnings(first(TRIAL_INDEX)),
      .groups = "drop"
    ) %>%
    group_by(!!trial_sym) %>%
    arrange(trial_index, desc(samples), Event) %>%
    slice(1) %>%
    pull(Event)
}

prepare_phase_samples <- function(data, cfg, participant_name, phase_label) {
  required_cols <- c("phase", "Event", "Time", "Pupil", cfg$video_col, cfg$trial_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols)) return(tibble())

  df <- data %>%
    filter(phase == cfg$phase_tag) %>%
    mutate(
      video_name = .data[[cfg$video_col]],
      trial_id = .data[[cfg$trial_col]]
    ) %>%
    filter(!is.na(video_name), video_name != ".", !is.na(trial_id))
  if (!nrow(df)) return(tibble())

  keep_events <- select_primary_events(df, cfg$trial_col)
  df <- df %>%
    filter(Event %in% keep_events) %>%
    arrange(Event, Time) %>%
    group_by(Event) %>%
    mutate(time_rel = Time - min(Time, na.rm = TRUE)) %>%
    ungroup()

  info_stim <- vapply(df$video_name, function(x) parse_video_info(x)$stimulus, character(1))
  info_suffix <- vapply(df$video_name, function(x) parse_video_info(x)$suffix, character(1))

  df <- df %>%
    mutate(
      stimulus = info_stim,
      suffix = info_suffix,
      participant = participant_name,
      phase_label = phase_label
    ) %>%
    filter(!is.na(stimulus), !is.na(suffix), suffix == cfg$suffix_keep)

  df
}

if (dir.exists(input_path)) {
  rds_files <- list.files(input_path, pattern = "^Data_.*\\.rds$", full.names = TRUE, ignore.case = TRUE)
  rds_files <- sort(rds_files)
} else if (file.exists(input_path)) {
  rds_files <- input_path
} else {
  stop("Input not found: ", input_path)
}
if (!length(rds_files)) stop("No RDS files found: ", input_path)

all_samples <- list()
for (rds_file in rds_files) {
  data <- readRDS(rds_file)
  if (!nrow(data)) next
  participant_name <- sub("^Data_", "", tools::file_path_sans_ext(basename(rds_file)))

  for (phase_label in names(phase_config)) {
    cfg <- phase_config[[phase_label]]
    phase_df <- prepare_phase_samples(data, cfg, participant_name, phase_label)
    if (!nrow(phase_df)) next
    all_samples[[length(all_samples) + 1]] <- phase_df
  }
}

if (!length(all_samples)) stop("No valid samples found.")
all_df <- bind_rows(all_samples)

all_df <- all_df %>%
  filter(is.finite(Pupil), is.finite(time_rel)) %>%
  mutate(
    time_bin = round(time_rel / bin_ms) * bin_ms
  ) %>%
  filter(time_bin >= 0, time_bin <= x_limit)

participant_bin <- all_df %>%
  group_by(phase_label, stimulus, participant, time_bin) %>%
  summarize(pupil_mean = mean(Pupil, na.rm = TRUE), .groups = "drop")

summary_df <- participant_bin %>%
  group_by(phase_label, stimulus, time_bin) %>%
  summarize(
    mean_pupil = mean(pupil_mean, na.rm = TRUE),
    sd_pupil = sd(pupil_mean, na.rm = TRUE),
    n = n_distinct(participant),
    .groups = "drop"
  ) %>%
  mutate(
    se_pupil = ifelse(n > 1, sd_pupil / sqrt(n), NA_real_),
    t_crit = ifelse(n > 1, stats::qt(0.975, df = n - 1), NA_real_),
    ci = se_pupil * t_crit,
    lower = mean_pupil - ci,
    upper = mean_pupil + ci
  )

line_colors <- c(white = "#FFD60A", black = "#0B3D91")

for (phase_label in names(phase_config)) {
  phase_df <- summary_df %>% filter(phase_label == !!phase_label)
  if (!nrow(phase_df)) next

  p <- ggplot(phase_df, aes(x = time_bin, y = mean_pupil, color = stimulus, fill = stimulus)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = c(7000, 25000), color = "black", linewidth = 0.3) +
    scale_color_manual(values = line_colors) +
    scale_fill_manual(values = line_colors) +
    scale_x_continuous(limits = c(0, x_limit), breaks = x_breaks, labels = function(x) format(x, big.mark = ",", trim = TRUE)) +
    labs(
      title = paste0(tools::toTitleCase(phase_label), " Phase (suffix ", phase_config[[phase_label]]$suffix_keep, ")"),
      x = "Time (ms)",
      y = "Pupil size (mm)",
      color = "Stimulus",
      fill = "Stimulus"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  if (phase_label %in% c("implicit", "explicit")) {
    p <- p + annotate(
      "rect",
      xmin = 7000,
      xmax = 12000,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey80",
      alpha = 0.25
    )
  }

  out_name <- paste0("pupil_mean_ci_", phase_label, ".png")
  ggsave(file.path(output_dir, out_name), p, width = 8, height = 4.5, dpi = 300)
  cat("✅ Saved:", file.path(output_dir, out_name), "\n")
}

cat("🎉 Mean pupil plots with CI saved to:", output_dir, "\n")
