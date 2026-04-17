library(dplyr)
library(ggplot2)
library(tidyr)

# Argument 1: RDS directory. Argument 2: output directory.
default_input <- "data/processed/RDS_Data"
default_output_dir <- "outputs/PupilSize/PupilMeanPhase"

time_window_ms <- c(7000, 12000)
sample_ms <- 40

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_dir <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

phase_config <- list(
  perception = list(phase_tag = "phase_p", video_col = "video_file_p", stage_keep = "0"),
  implicit   = list(phase_tag = "phase_i", video_col = "video_file_i", stage_keep = "1"),
  explicit   = list(phase_tag = "phase_e", video_col = "video_file_e", stage_keep = "1")
)

parse_video_meta <- function(video_name) {
  base <- tools::file_path_sans_ext(basename(video_name))
  pieces <- strsplit(base, "-")[[1]]
  color <- if (length(pieces)) tolower(pieces[1]) else NA_character_
  stage <- if (length(pieces) >= 2) pieces[2] else "0"
  list(color = color, stage = stage)
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

# Collect participant-level pupil traces for each phase and stimulus combination.
all_rows <- list()

for (rds_file in rds_files) {
  df <- readRDS(rds_file)
  if (!nrow(df)) next

  participant_id <- NA_character_
  if ("Subject" %in% names(df)) {
    participant_id <- pad_id(na.omit(as.character(df$Subject))[1])
  }
  if (is.na(participant_id) || !nzchar(participant_id)) {
    participant_id <- pad_id(sub("^Data_", "", tools::file_path_sans_ext(basename(rds_file))))
  }

  for (phase_label in names(phase_config)) {
    cfg <- phase_config[[phase_label]]
    if (!(cfg$phase_tag %in% df$phase)) next
    if (!(cfg$video_col %in% names(df))) next

    phase_df <- df %>%
      filter(.data$phase == cfg$phase_tag, !is.na(Pupil), !is.na(Time), !is.na(Event)) %>%
      mutate(video_file = .data[[cfg$video_col]]) %>%
      filter(!is.na(video_file), video_file != ".")

    if (!nrow(phase_df)) next

    meta <- lapply(phase_df$video_file, parse_video_meta)
    phase_df$stimulus <- vapply(meta, function(x) x$color, character(1))
    phase_df$stage <- vapply(meta, function(x) x$stage, character(1))

    phase_df <- phase_df %>%
      filter(!is.na(stimulus), stage == cfg$stage_keep) %>%
      arrange(Event, Time) %>%
      group_by(Event) %>%
      mutate(time_rel = Time - min(Time, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(time_bin = round(time_rel / sample_ms) * sample_ms)

    if (!nrow(phase_df)) next

    part_df <- phase_df %>%
      group_by(participant = participant_id, phase = phase_label, stimulus, time_bin) %>%
      summarize(pupil = mean(Pupil, na.rm = TRUE), .groups = "drop")

    all_rows[[length(all_rows) + 1]] <- part_df
  }
}

# Average across participants and compute confidence intervals for each time bin.
all_df <- bind_rows(all_rows)
if (!nrow(all_df)) stop("No pupil data for plots.")

summary_df <- all_df %>%
  group_by(phase, stimulus, time_bin) %>%
  summarize(
    mean = mean(pupil, na.rm = TRUE),
    sd = sd(pupil, na.rm = TRUE),
    n = n_distinct(participant),
    .groups = "drop"
  ) %>%
  mutate(
    se = ifelse(n > 1, sd / sqrt(n), NA_real_),
    t_crit = ifelse(n > 1, stats::qt(0.975, df = n - 1), NA_real_),
    ci = se * t_crit
  ) %>%
  mutate(
    stimulus = factor(stimulus, levels = c("white", "black"))
  )

# Plot one phase-specific mean pupil trace with confidence intervals.
plot_phase <- function(phase_label, add_window, output_file) {
  df_phase <- summary_df %>%
    filter(phase == phase_label)
  if (!nrow(df_phase)) return()

  phase_titles <- c(
    perception = "Perception",
    explicit = "Instructed Imagery",
    implicit = "Passive Viewing"
  )
  phase_title <- phase_titles[[phase_label]]

  p <- ggplot(df_phase, aes(x = time_bin, y = mean, color = stimulus, fill = stimulus)) +
    {if (add_window) annotate(
      "rect",
      xmin = time_window_ms[1],
      xmax = time_window_ms[2],
      ymin = -Inf,
      ymax = Inf,
      fill = "grey80",
      alpha = 0.25
    )} +
    geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.15, color = NA) +
    geom_vline(xintercept = c(7000, 25000), color = "black", linewidth = 0.4) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(white = "#FFD60A", black = "#0B3D91")) +
    scale_fill_manual(values = c(white = "#FFD60A", black = "#0B3D91")) +
    scale_x_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000)) +
    labs(
      x = "Time (ms)",
      y = "Pupil (mean)",
      color = "Stimulus",
      fill = "Stimulus",
      caption = phase_title
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      plot.caption = element_text(hjust = 0.5, face = "bold"),
      plot.caption.position = "plot",
      text = element_text(family = "Times New Roman")
    )

  ggsave(output_file, p, width = 7, height = 4.5, dpi = 300)
  cat("Saved:", output_file, "\n")
}

# Export one figure per phase using the appropriate occlusion window setting.
plot_phase("perception", FALSE, file.path(output_dir, "pupil_mean_perception.png"))
plot_phase("explicit", TRUE, file.path(output_dir, "pupil_mean_explicit.png"))
plot_phase("implicit", TRUE, file.path(output_dir, "pupil_mean_implicit.png"))
