library(dplyr)
library(readxl)
library(stringr)
library(tidyr)

# === Config ===
default_rds_dir <- "data/processed/RDS_Data"
output_csv <- "data/processed/gaze_target_cosine_occlusion.csv"

# Screen size for plotting (pixels)
screen_width <- 4480
screen_height <- 2520

# Coin path (Unity export)
coin_path_white_file <- "data/reference/coin_path_white.xlsx"
coin_path_black_file <- "data/reference/coin_path_black.xlsx"
coin_path_screen_w <- screen_width
coin_path_screen_h <- screen_height

# Do not swap coin paths
swap_coin_paths <- FALSE

# Use a fixed time window (ms)
time_window_ms <- c(7000, 12000)

# === Helpers ===
pad_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.0+$", "", x)
  suppressWarnings({
    x_num <- as.integer(x)
  })
  ifelse(!is.na(x_num), sprintf("%08d", x_num), x)
}

parse_video_meta <- function(video_file) {
  base <- tools::file_path_sans_ext(video_file)
  pieces <- str_split(base, "-", simplify = TRUE)
  color <- if (ncol(pieces) >= 1) tolower(pieces[1]) else "unknown"
  stage <- if (ncol(pieces) >= 2) pieces[2] else "0"
  list(color = color, stage = stage)
}

load_coin_path <- function(path, screen_width, screen_height, base_w, base_h) {
  if (!file.exists(path)) return(NULL)
  df <- readxl::read_excel(path)
  if (!all(c("screen_x", "screen_y") %in% names(df))) {
    warning("Coin path missing screen_x/screen_y: ", path)
    return(NULL)
  }
  if ("t" %in% names(df)) {
    df <- df %>% arrange(t)
  }
  base_w <- as.numeric(base_w)
  base_h <- as.numeric(base_h)
  if (is.na(base_w) || base_w <= 0) base_w <- screen_width
  if (is.na(base_h) || base_h <= 0) base_h <- screen_height
  sx <- screen_width / base_w
  sy <- screen_height / base_h
  df %>%
    transmute(
      t = as.numeric(t),
      screen_x = as.numeric(screen_x) * sx,
      screen_y = as.numeric(screen_y) * sy
    ) %>%
    drop_na(t, screen_x, screen_y)
}

interp_target <- function(time_s, coin_df) {
  if (is.null(coin_df) || !nrow(coin_df)) {
    return(list(x = rep(NA_real_, length(time_s)), y = rep(NA_real_, length(time_s))))
  }
  x <- approx(coin_df$t, coin_df$screen_x, xout = time_s, rule = 2)$y
  y <- approx(coin_df$t, coin_df$screen_y, xout = time_s, rule = 2)$y
  list(x = x, y = y)
}

safe_cosine <- function(vgx, vgy, vtx, vty) {
  denom <- sqrt(vgx^2 + vgy^2) * sqrt(vtx^2 + vty^2)
  cos_theta <- (vgx * vtx + vgy * vty) / denom
  cos_theta[is.nan(cos_theta) | is.infinite(cos_theta)] <- NA_real_
  cos_theta
}

# === Load Assets ===
# Load the target motion paths used to define the expected movement direction.
coin_paths <- list(
  white = load_coin_path(
    if (swap_coin_paths) coin_path_black_file else coin_path_white_file,
    screen_width, screen_height, coin_path_screen_w, coin_path_screen_h
  ),
  black = load_coin_path(
    if (swap_coin_paths) coin_path_white_file else coin_path_black_file,
    screen_width, screen_height, coin_path_screen_w, coin_path_screen_h
  )
)

# === Inputs ===
args <- commandArgs(trailingOnly = TRUE)
rds_dir <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_rds_dir

rds_files <- list.files(rds_dir, pattern = "^Data_.*\\.rds$", full.names = TRUE)
if (!length(rds_files)) stop("No RDS files found: ", rds_dir)

phase_configs <- list(
  perception = list(phase_tag = "phase_p", video_col = "video_file_p", trial_col = "trial_p"),
  implicit   = list(phase_tag = "phase_i", video_col = "video_file_i", trial_col = "trial_i"),
  explicit   = list(phase_tag = "phase_e", video_col = "video_file_e", trial_col = "trial_e")
)

phase_map <- c(
  phase_p = "perception",
  phase_i = "implicit",
  phase_e = "explicit"
)

# Iterate over participant RDS files and compute trial-level gaze-target cosine alignment.
all_results <- list()

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

  for (cfg in phase_configs) {
    phase_tag <- cfg$phase_tag
    video_col <- cfg$video_col
    trial_col <- cfg$trial_col
    if (!(phase_tag %in% names(phase_map))) next
    if (!(phase_tag %in% df$phase)) next
    if (!(video_col %in% names(df))) next
    if (!(trial_col %in% names(df))) next

    phase_df <- df %>%
      filter(.data$phase == phase_tag, !is.na(Time), !is.na(Gaze_X), !is.na(Gaze_Y)) %>%
      mutate(
        trial_id = .data[[trial_col]],
        video_file = .data[[video_col]]
      ) %>%
      filter(!is.na(video_file))

    if (!nrow(phase_df)) next

    trial_groups <- phase_df %>% group_by(Event) %>% group_split()

    for (trial_df in trial_groups) {
      video_file <- na.omit(trial_df$video_file)[1]
      if (is.na(video_file) || !nzchar(video_file)) next
      meta <- parse_video_meta(video_file)
      if (phase_tag %in% c("phase_i", "phase_e") && meta$stage != "1") next

      coin_df <- coin_paths[[meta$color]]

      trial_df <- trial_df %>%
        arrange(Time) %>%
        mutate(
          time_s = Time / 1000
        )

      target_pos <- interp_target(trial_df$time_s, coin_df)
      trial_df$target_x <- target_pos$x
      trial_df$target_y <- target_pos$y

      # Convert gaze and target motion into velocity vectors before computing cosine similarity.
      trial_df <- trial_df %>%
        mutate(
          dt_s = (Time - lag(Time)) / 1000,
          vgx = (Gaze_X - lag(Gaze_X)) / dt_s,
          vgy = (Gaze_Y - lag(Gaze_Y)) / dt_s,
          vtx = (target_x - lag(target_x)) / dt_s,
          vty = (target_y - lag(target_y)) / dt_s,
          cos_theta = safe_cosine(vgx, vgy, vtx, vty)
        )

      trial_df <- trial_df %>%
        mutate(
          in_window = Time >= time_window_ms[1] & Time <= time_window_ms[2]
        )

      # Summarize cosine alignment within the occlusion window for each trial.
      window_df <- trial_df %>% filter(in_window)
      if (!nrow(window_df)) {
        result <- data.frame(
          participant = participant_id,
          phase = phase_map[[phase_tag]],
          stimulus = meta$color,
          stage = meta$stage,
          event = as.character(trial_df$Event[1]),
          trial_id = trial_df$trial_id[1],
          cos_mean = NA_real_,
          cos_median = NA_real_,
          n_samples = 0,
          occlusion_start_ms = NA_real_,
          occlusion_end_ms = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        result <- data.frame(
          participant = participant_id,
          phase = phase_map[[phase_tag]],
          stimulus = meta$color,
          stage = meta$stage,
          event = as.character(trial_df$Event[1]),
          trial_id = trial_df$trial_id[1],
          cos_mean = mean(window_df$cos_theta, na.rm = TRUE),
          cos_median = median(window_df$cos_theta, na.rm = TRUE),
          n_samples = sum(!is.na(window_df$cos_theta)),
          occlusion_start_ms = min(window_df$Time, na.rm = TRUE),
          occlusion_end_ms = max(window_df$Time, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
      all_results[[length(all_results) + 1]] <- result
    }
  }
}

# Save the combined trial-level cosine summary table.
results_df <- bind_rows(all_results)
if (!nrow(results_df)) stop("No occlusion cosine results generated. Check inputs/paths.")

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(results_df, output_csv, row.names = FALSE)
cat("Saved:", output_csv, "\n")
