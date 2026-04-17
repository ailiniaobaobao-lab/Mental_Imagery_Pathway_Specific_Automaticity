library(dplyr) 
library(tidyr) 
library(lme4) 
library(readxl) 
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest) 
}

# input
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_dorsal <- "data/processed/dorsal_trial_errors.csv"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE) 
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
dorsal_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_dorsal
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE) 

output_txt <- file.path(output_dir, "automaticity_mixed_models.txt") 

strength_eps <- 1e-6 
rmse_eps <- 1e-6 


pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

exclude_participants <- pad_id(exclude_participants)


participants_info <- readxl::read_excel("data/metadata/Participants_info.xlsx")
participants_info <- participants_info %>%
  transmute(
    participant = pad_id(`subject#`),
    stim_only = as.character(`...10`)
  )
black_only_ids <- participants_info %>%
  filter(tolower(stim_only) == "black only") %>%
  pull(participant)
white_only_ids <- participants_info %>%
  filter(tolower(stim_only) == "white only") %>%
  pull(participant)

apply_stimulus_restrictions <- function(df) {
  df %>%
    filter(!(participant %in% black_only_ids & stimulus == "white")) %>%
    filter(!(participant %in% white_only_ids & stimulus == "black"))
}

# ==== Load Data ====
ventral <- read.csv(ventral_path, stringsAsFactors = FALSE)
dorsal <- read.csv(dorsal_path, stringsAsFactors = FALSE)

# ==== Ventral: trial-level mixed model ====
ventral_df_all <- ventral %>%
  mutate(
    phase = case_when(
      phase == "phase_i" ~ "implicit",
      phase == "phase_e" ~ "explicit",
      phase == "phase_p" ~ "perception",
      TRUE ~ NA_character_
    ),
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    stimulus = factor(stimulus, levels = c("white", "black")),
    participant = pad_id(participant)
  ) %>%
  filter(!is.na(pupil_response), !is.na(participant), !is.na(stimulus), !is.na(phase)) %>%
  filter(!participant %in% exclude_participants)

ventral_df_model <- ventral_df_all %>%
  filter(!participant %in% black_only_ids, !participant %in% white_only_ids) %>%
  mutate(participant = factor(participant))

ventral_perc_stats_z_raw <- ventral_df_model %>%
  filter(phase == "perception") %>%
  group_by(participant, stimulus) %>%
  summarize(
    base_mean = mean(pupil_response, na.rm = TRUE),
    base_sd = sd(pupil_response, na.rm = TRUE),
    .groups = "drop"
  )

ventral_z_all_raw <- ventral_df_model %>%
  left_join(ventral_perc_stats_z_raw, by = c("participant", "stimulus")) %>%
  mutate(dv_z = (pupil_response - base_mean) / base_sd) %>%
  filter(is.finite(dv_z))

ventral_z_df_raw <- ventral_z_all_raw %>%
  filter(phase %in% c("explicit", "implicit"))

ventral_model_z <- lmer(dv_z ~ phase * stimulus + (1 | participant), data = ventral_z_df_raw, REML = FALSE)

# ==== Dorsal: trial-level mixed model (log RMSE) ====
dorsal_df_all <- dorsal %>%
  mutate(
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    stimulus = factor(stimulus, levels = c("white", "black")),
    traj_rmse_log = log(pmax(as.numeric(traj_rmse), rmse_eps)),
    participant = pad_id(participant)
  ) %>%
  filter(is.finite(traj_rmse_log), !is.na(participant), !is.na(phase), !is.na(stimulus)) %>%
  filter(!participant %in% exclude_participants)

dorsal_df_model <- dorsal_df_all %>%
  mutate(participant = factor(participant))

# ==== Combined: trial-level (log RMSE) ====
dorsal_rmse_all <- dorsal_df_all %>%
  filter(is.finite(traj_rmse_log), !is.na(participant), !is.na(phase), !is.na(stimulus))

ventral_signed_all <- ventral_df_all %>%
  transmute(
    participant,
    phase,
    stimulus,
    pathway = "ventral",
    dv_raw = ifelse(stimulus == "white", -pupil_response, pupil_response)
  )

ventral_perc_stats_z_signed <- ventral_signed_all %>%
  filter(phase == "perception") %>%
  group_by(participant, stimulus) %>%
  summarize(
    base_mean = mean(dv_raw, na.rm = TRUE),
    base_sd = sd(dv_raw, na.rm = TRUE),
    .groups = "drop"
  )

ventral_z_all_signed <- ventral_signed_all %>%
  left_join(ventral_perc_stats_z_signed, by = c("participant", "stimulus")) %>%
  mutate(dv_z = (dv_raw - base_mean) / base_sd) %>%
  filter(is.finite(dv_z))

dorsal_perc_stats_z <- dorsal_rmse_all %>%
  filter(phase == "perception") %>%
  group_by(participant, stimulus) %>%
  summarize(
    base_mean = mean(traj_rmse_log, na.rm = TRUE),
    base_sd = sd(traj_rmse_log, na.rm = TRUE),
    .groups = "drop"
  )

dorsal_z_all <- dorsal_rmse_all %>%
  left_join(dorsal_perc_stats_z, by = c("participant", "stimulus")) %>%
  mutate(dv_z = (traj_rmse_log - base_mean) / base_sd) %>%
  filter(is.finite(dv_z))

dorsal_z_df <- dorsal_z_all %>%
  filter(phase %in% c("explicit", "implicit"))

dorsal_model_z <- lmer(dv_z ~ phase + (1 | participant), data = dorsal_z_df %>% mutate(participant = factor(participant)), REML = FALSE)

build_combined <- function(ventral_df, dorsal_df, phases) {
  ventral_use <- ventral_df %>%
    filter(phase %in% phases) %>%
    transmute(
      participant,
      phase,
      stimulus,
      pathway = "ventral",
      dv_z
    )
  dorsal_use <- dorsal_df %>%
    filter(phase %in% phases) %>%
    transmute(
      participant,
      phase,
      stimulus,
      pathway = "dorsal",
      dv_z
    )
  combined <- bind_rows(ventral_use, dorsal_use) %>%
    mutate(
      participant = factor(participant),
      pathway = factor(pathway, levels = c("ventral", "dorsal")),
      phase = factor(phase, levels = phases),
      stimulus = factor(stimulus, levels = c("white", "black"))
    )
  combined
}

ventral_z_all_signed <- apply_stimulus_restrictions(ventral_z_all_signed)
dorsal_z_all <- apply_stimulus_restrictions(dorsal_z_all)

combined_signed_rmse <- build_combined(ventral_z_all_signed, dorsal_z_all, c("explicit", "implicit"))
combined_model_signed_rmse <- lmer(dv_z ~ phase * pathway + stimulus + (1 | participant), data = combined_signed_rmse, REML = FALSE)

# ==== Dorsal time-shift control (if available) ====
shift_cols <- grep("^traj_rmse_shift", names(dorsal_df_all), value = TRUE)
shift_control_text <- character()
if (length(shift_cols) >= 1) {
  shift_col <- shift_cols[1]
  dorsal_shift_df <- dorsal_df_all %>%
    filter(is.finite(.data[[shift_col]])) %>%
    transmute(
      participant = factor(participant),
      phase,
      stimulus,
      condition = "real",
      rmse = as.numeric(traj_rmse)
    )
  dorsal_shift_df <- bind_rows(
    dorsal_shift_df,
    dorsal_df_all %>%
      filter(is.finite(.data[[shift_col]])) %>%
      transmute(
        participant = factor(participant),
        phase,
        stimulus,
        condition = "shift",
        rmse = as.numeric(.data[[shift_col]])
      )
  ) %>%
    filter(is.finite(rmse)) %>%
    mutate(
      condition = factor(condition, levels = c("real", "shift")),
      phase = factor(phase, levels = c("perception", "explicit", "implicit")),
      stimulus = factor(stimulus, levels = c("white", "black")),
      tracking_strength = -log(pmax(rmse, rmse_eps))
    )

  shift_means <- dorsal_shift_df %>%
    group_by(phase, condition) %>%
    summarize(mean_strength = mean(tracking_strength, na.rm = TRUE), .groups = "drop") %>%
    arrange(phase, condition)

  shift_model <- lmer(tracking_strength ~ condition * phase + stimulus + (1 | participant),
                      data = dorsal_shift_df, REML = FALSE)

  shift_control_text <- c(
    "",
    "=== Dorsal time-shift control (target delayed 2000ms; tracking_strength = -log(RMSE)) ===",
    "Per-phase mean tracking_strength (higher = better):",
    capture.output(shift_means),
    "",
    "Mixed model: tracking_strength ~ condition * phase + stimulus + (1|participant)",
    capture.output(summary(shift_model))
  )
}

# output
summary_text <- c(
  "=== Ventral: pupil_response (z-scored vs perception) ~ phase * stimulus + (1|participant) ===",
  capture.output(summary(ventral_model_z)),
  "",
  "=== Dorsal: log(RMSE) (z-scored vs perception) ~ phase + (1|participant) ===",
  capture.output(summary(dorsal_model_z)),
  "",
  "=== Combined signed (explicit/implicit): dv_z ~ phase * pathway + stimulus + (1|participant) [log(RMSE) + perception-z; white flipped] ===",
  capture.output(summary(combined_model_signed_rmse)),
  shift_control_text
)

writeLines(summary_text, con = output_txt)

cat("save：", output_txt, "\n")
