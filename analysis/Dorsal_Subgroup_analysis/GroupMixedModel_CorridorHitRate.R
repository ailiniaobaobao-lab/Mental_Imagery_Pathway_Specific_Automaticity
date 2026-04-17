library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(lme4)
library(ggplot2)
library(emmeans)
library(eyetrackingR)  
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}

# input
default_rds_dir <- "data/processed/RDS_Data"
default_info <- "data/metadata/Participants_info.xlsx"
default_output_dir <- "outputs/Gaze_corrodor"

# Coin path (Unity export)
coin_path_white_file <- "data/reference/coin_path_white.xlsx"
coin_path_black_file <- "data/reference/coin_path_black.xlsx"

# Screen size for plotting (pixels)
screen_width <- 4480
screen_height <- 2520

# Use a fixed time window (ms)
time_window_ms <- c(7000, 12000)

# Time-shift control (ms): target trajectory delayed by 2s
time_shift_ms <- 2000

# Do not swap coin paths
swap_coin_paths <- FALSE

# Flip coin Y to match gaze screen coordinates
flip_coin_y <- TRUE

args <- commandArgs(trailingOnly = TRUE)
rds_dir <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_rds_dir
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

parse_video_meta <- function(video_file) {
  base <- tools::file_path_sans_ext(video_file)
  pieces <- str_split(base, "-", simplify = TRUE)
  color <- if (ncol(pieces) >= 1) tolower(pieces[1]) else "unknown"
  stage <- if (ncol(pieces) >= 2) pieces[2] else "0"
  list(color = color, stage = stage)
}

load_coin_path <- function(path, screen_width, screen_height, flip_y) {
  if (!file.exists(path)) return(NULL)
  df <- readxl::read_excel(path)
  if (!all(c("screen_x", "screen_y") %in% names(df))) {
    warning("Coin path missing screen_x/screen_y: ", path)
    return(NULL)
  }
  if ("t" %in% names(df)) {
    df <- df %>% arrange(t)
  }
  out <- df %>%
    transmute(
      t = as.numeric(t),
      screen_x = as.numeric(screen_x),
      screen_y = as.numeric(screen_y)
    ) %>%
    drop_na(t, screen_x, screen_y)
  if (isTRUE(flip_y)) {
    out <- out %>% mutate(screen_y = screen_height - screen_y)
  }
  out
}

interp_target <- function(time_s, coin_df) {
  if (is.null(coin_df) || !nrow(coin_df)) {
    return(list(x = rep(NA_real_, length(time_s)), y = rep(NA_real_, length(time_s))))
  }
  x <- approx(coin_df$t, coin_df$screen_x, xout = time_s, rule = 2)$y
  y <- approx(coin_df$t, coin_df$screen_y, xout = time_s, rule = 2)$y
  list(x = x, y = y)
}

coin_paths <- list(
  white = load_coin_path(
    if (swap_coin_paths) coin_path_black_file else coin_path_white_file,
    screen_width, screen_height, flip_coin_y
  ),
  black = load_coin_path(
    if (swap_coin_paths) coin_path_white_file else coin_path_black_file,
    screen_width, screen_height, flip_coin_y
  )
)

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

rds_files <- list.files(rds_dir, pattern = "^Data_.*\\.rds$", full.names = TRUE)
if (!length(rds_files)) stop("No RDS files found: ", rds_dir)

sample_rows <- list()
trial_rows <- list()

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
      if (is.null(coin_df) || !nrow(coin_df)) next
      
      trial_df <- trial_df %>%
        arrange(Time) %>%
        mutate(time_s = Time / 1000)
      
      target_pos <- interp_target(trial_df$time_s, coin_df)
      target_pos_shift <- interp_target(trial_df$time_s - time_shift_ms / 1000, coin_df)
      trial_df$target_x <- target_pos$x
      trial_df$target_y <- target_pos$y
      trial_df$target_x_shift <- target_pos_shift$x
      trial_df$target_y_shift <- target_pos_shift$y
      
      trial_df <- trial_df %>%
        mutate(
          dist = sqrt((Gaze_X - target_x)^2 + (Gaze_Y - target_y)^2),
          dist_shift = sqrt((Gaze_X - target_x_shift)^2 + (Gaze_Y - target_y_shift)^2),
          in_window = Time >= time_window_ms[1] & Time <= time_window_ms[2]
        )
      
      window_df <- trial_df %>% filter(in_window)
      
      trial_key <- data.frame(
        participant = participant_id,
        phase = phase_map[[phase_tag]],
        stimulus = meta$color,
        stage = meta$stage,
        event = as.character(trial_df$Event[1]),
        trial_id = trial_df$trial_id[1],
        stringsAsFactors = FALSE
      )
      trial_rows[[length(trial_rows) + 1]] <- trial_key
      
      if (nrow(window_df)) {
        sample_rows[[length(sample_rows) + 1]] <- window_df %>%
          transmute(
            participant = participant_id,
            phase = phase_map[[phase_tag]],
            stimulus = meta$color,
            stage = meta$stage,
            event = as.character(Event),
            trial_id = trial_id,
            Time = Time,
            gaze_x = Gaze_X,
            gaze_y = Gaze_Y,
            target_x = target_x,
            target_y = target_y,
            dist = dist,
            condition = "real"
          )
        sample_rows[[length(sample_rows) + 1]] <- window_df %>%
          transmute(
            participant = participant_id,
            phase = phase_map[[phase_tag]],
            stimulus = meta$color,
            stage = meta$stage,
            event = as.character(Event),
            trial_id = trial_id,
            Time = Time,
            gaze_x = Gaze_X,
            gaze_y = Gaze_Y,
            target_x = target_x_shift,
            target_y = target_y_shift,
            dist = dist_shift,
            condition = "shift"
          )
      }
    }
  }
}

samples_all <- bind_rows(sample_rows)
trial_catalog <- bind_rows(trial_rows) %>%
  distinct(participant, phase, stimulus, stage, event, trial_id)

if (!nrow(samples_all)) stop("No samples found in the time window.")

samples_df <- samples_all %>% filter(condition == "real")


# corridor：perception 95% = R

thresholds <- samples_df %>%
  filter(phase == "perception") %>%
  group_by(participant, stimulus) %>%
  summarize(
    threshold_95 = quantile(dist, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# dynamic AOI bbox
aoi_df <- samples_all %>%
  left_join(thresholds, by = c("participant", "stimulus")) %>%
  mutate(
    x_min = target_x - threshold_95,
    x_max = target_x + threshold_95,
    y_min = target_y - threshold_95,
    y_max = target_y + threshold_95
  ) %>%

  select(participant, phase, stimulus, stage, event, trial_id, Time, condition,
         x_min, x_max, y_min, y_max)

hit_samples <- eyetrackingR::add_aoi(
  data = samples_all,
  aoi_dataframe = aoi_df,
  x_col = "gaze_x",
  y_col = "gaze_y",
  aoi_name = "Corridor",
  x_min_col = "x_min",
  x_max_col = "x_max",
  y_min_col = "y_min",
  y_max_col = "y_max"
)

# trial-level hit rate
trial_hits <- hit_samples %>%
  group_by(participant, phase, stimulus, stage, event, trial_id) %>%
  summarize(
    hit_rate = mean(Corridor, na.rm = TRUE),
    n_samples = sum(!is.na(Corridor)),
    .groups = "drop"
  )

trial_out <- trial_catalog %>%
  left_join(trial_hits, by = c("participant", "phase", "stimulus", "stage", "event", "trial_id")) %>%
  left_join(thresholds, by = c("participant", "stimulus"))

out_trials <- file.path(output_dir, "gaze_corridor_hit_rate_trials.csv")
write.csv(trial_out, out_trials, row.names = FALSE)

# time-shift control output (real vs shift)
trial_hits_all <- hit_samples %>%
  group_by(participant, phase, stimulus, stage, event, trial_id, condition) %>%
  summarize(
    hit_rate = mean(Corridor, na.rm = TRUE),
    n_samples = sum(!is.na(Corridor)),
    .groups = "drop"
  )

trial_catalog_all <- tidyr::expand_grid(
  trial_catalog,
  condition = c("real", "shift")
)

trial_out_all <- trial_catalog_all %>%
  left_join(trial_hits_all, by = c("participant", "phase", "stimulus", "stage", "event", "trial_id", "condition")) %>%
  left_join(thresholds, by = c("participant", "stimulus"))

out_trials_shift <- file.path(output_dir, "gaze_corridor_hit_rate_trials_timeshift.csv")
write.csv(trial_out_all, out_trials_shift, row.names = FALSE)

shift_summary <- trial_out_all %>%
  filter(is.finite(hit_rate)) %>%
  group_by(phase, condition) %>%
  summarize(mean_hit_rate = mean(hit_rate, na.rm = TRUE), .groups = "drop") %>%
  arrange(phase, condition)

shift_model_df <- trial_out_all %>%
  filter(is.finite(hit_rate)) %>%
  mutate(
    participant = factor(participant),
    phase = factor(phase, levels = c("perception", "explicit", "implicit")),
    stimulus = factor(stimulus, levels = c("white", "black")),
    condition = factor(condition, levels = c("real", "shift"))
  )

shift_model <- lmer(hit_rate ~ condition * phase + stimulus + (1 | participant), data = shift_model_df, REML = FALSE)

out_txt_shift <- file.path(output_dir, "gaze_corridor_hit_rate_timeshift.txt")
shift_text <- c(
  "=== Corridor hit rate time-shift control (target delayed 2000ms) ===",
  "Per-phase mean hit_rate (higher = better):",
  capture.output(shift_summary),
  "",
  "Mixed model: hit_rate ~ condition * phase + stimulus + (1|participant)",
  capture.output(summary(shift_model))
)
writeLines(shift_text, con = out_txt_shift)

# Explicit/Implicit only (no perception)
shift_model_df_epi <- shift_model_df %>%
  filter(phase %in% c("explicit", "implicit")) %>%
  droplevels()

shift_summary_epi <- shift_model_df_epi %>%
  group_by(phase, condition) %>%
  summarize(mean_hit_rate = mean(hit_rate, na.rm = TRUE), .groups = "drop") %>%
  arrange(phase, condition)

shift_model_epi <- lmer(hit_rate ~ condition * phase + stimulus + (1 | participant),
                        data = shift_model_df_epi, REML = FALSE)

out_txt_shift_epi <- file.path(output_dir, "gaze_corridor_hit_rate_timeshift_explicit_implicit.txt")
shift_text_epi <- c(
  "=== Corridor hit rate time-shift control (explicit/implicit only) ===",
  "Per-phase mean hit_rate (higher = better):",
  capture.output(shift_summary_epi),
  "",
  "Mixed model: hit_rate ~ condition * phase + stimulus + (1|participant)",
  capture.output(summary(shift_model_epi))
)
writeLines(shift_text_epi, con = out_txt_shift_epi)

# group map
info <- readxl::read_excel(info_path)
col_h <- names(info)[8]
group_map <- info %>%
  transmute(
    participant = pad_id(.data[["subject#"]]),
    group = as.character(.data[[col_h]])
  ) %>%
  filter(group %in% c("B", "N", "E"))

model_df <- trial_out %>%
  inner_join(group_map, by = "participant") %>%
  filter(is.finite(hit_rate)) %>%
  mutate(
    participant = factor(participant),
    phase = factor(phase, levels = c("perception", "explicit", "implicit")),
    group = factor(group, levels = c("B", "N", "E")),
    stimulus = factor(stimulus, levels = c("white", "black"))
  )

model <- lmer(hit_rate ~ phase * group + stimulus + (1 | participant), data = model_df, REML = FALSE)

out_txt <- file.path(output_dir, "gaze_corridor_hit_rate_mixed_model.txt")
summary_text <- c(
  "=== Corridor hit rate mixed model (AOI): hit_rate ~ phase * group + stimulus + (1|participant) ===",
  capture.output(summary(model))
)
writeLines(summary_text, con = out_txt)

emm_df <- as.data.frame(emmeans(model, ~ phase * group)) %>%
  mutate(section = "emmeans_phase_group")
out_emm <- file.path(output_dir, "gaze_corridor_hit_rate_emmeans.csv")
write.csv(emm_df, out_emm, row.names = FALSE)

plot_df <- emm_df %>%
  mutate(
    phase = factor(phase, levels = c("perception", "explicit", "implicit")),
    group = factor(group, levels = c("B", "E", "N"))
  )

phase_labels <- c(
  perception = "Perception",
  explicit = "Instructed Imagery",
  implicit = "Passive Viewing"
)

group_labels <- c(
  B = "B",
  E = "I",
  N = "N"
)

p <- ggplot(plot_df, aes(x = phase, y = emmean, color = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.08, linewidth = 0.6) +
  labs(
    x = "Phase",
    y = "Hit rate",
    color = "Group"
  ) +
  scale_x_discrete(labels = phase_labels) +
  scale_color_discrete(labels = group_labels) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "Times New Roman"))

out_plot <- file.path(output_dir, "gaze_corridor_hit_rate_emmeans.png")
ggsave(out_plot, p, width = 6.5, height = 4.5, dpi = 300)

cat("Saved:", out_trials, "\n")
cat("Saved:", out_txt, "\n")
cat("Saved:", out_emm, "\n")
cat("Saved:", out_plot, "\n")
