library(dplyr)
library(tidyr)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}

# Argument 1: ventral trial CSV. Argument 2: corridor trial CSV. Argument 3: optional output directory.
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_corridor <- "data/processed/gaze_corridor_hit_rate_trials.csv"
default_output_dir <- "outputs/Gaze_corrodor"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
corridor_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_corridor
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_txt <- file.path(output_dir, "automaticity_mixed_models_corridor_unified.txt")

exclude_participants <- c("11191401", "11191503", "11210905", "11241311", "11241210", "03041338")

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

exclude_participants <- pad_id(exclude_participants)

# ==== Load Data ====
ventral <- read.csv(ventral_path, stringsAsFactors = FALSE)
corridor <- read.csv(corridor_path, stringsAsFactors = FALSE)

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

# ==== Corridor (dorsal replacement): trial-level mixed model (raw hit_rate) ====
corridor_df_all <- corridor %>%
  mutate(
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    stimulus = factor(stimulus, levels = c("white", "black")),
    participant = pad_id(participant)
  ) %>%
  filter(!is.na(hit_rate), !is.na(participant), !is.na(phase), !is.na(stimulus)) %>%
  filter(!participant %in% exclude_participants)

corridor_raw_df <- corridor_df_all %>%
  filter(phase %in% c("explicit", "implicit"))

corridor_model_raw <- lmer(hit_rate ~ phase + stimulus + (1 | participant), data = corridor_raw_df %>% mutate(participant = factor(participant)), REML = FALSE)

# ==== Combined: ventral signed + corridor raw (unified scale within participant) ====
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

corridor_for_combined <- corridor_df_all %>%
  transmute(
    participant,
    phase,
    stimulus,
    pathway = "dorsal",
    dv_raw = hit_rate
  )

# Build the shared long-format dataset used for the cross-pathway mixed model.
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
      dv_raw
    )
  combined <- bind_rows(ventral_use %>% rename(dv_raw = dv_z), dorsal_use) %>%
    group_by(participant) %>%
    mutate(
      dv_mean = mean(dv_raw, na.rm = TRUE),
      dv_sd = sd(dv_raw, na.rm = TRUE),
      dv_z = (dv_raw - dv_mean) / dv_sd
    ) %>%
    ungroup() %>%
    filter(is.finite(dv_z)) %>%
    mutate(
      participant = factor(participant),
      pathway = factor(pathway, levels = c("ventral", "dorsal")),
      phase = factor(phase, levels = phases),
      stimulus = factor(stimulus, levels = c("white", "black"))
    )
  combined
}

combined_signed_corridor <- build_combined(ventral_z_all_signed, corridor_for_combined, c("explicit", "implicit"))
combined_model_signed_corridor <- lmer(dv_z ~ phase * pathway + stimulus + (1 | participant), data = combined_signed_corridor, REML = FALSE)

# ==== Output ====
summary_text <- c(
  "=== Ventral: pupil_response (z-scored vs perception) ~ phase * stimulus + (1|participant) ===",
  capture.output(summary(ventral_model_z)),
  "",
  "=== Corridor (dorsal replacement): hit_rate (raw) ~ phase + stimulus + (1|participant) ===",
  capture.output(summary(corridor_model_raw)),
  "",
  "=== Combined signed (explicit/implicit): dv_z ~ phase * pathway + stimulus + (1|participant) [participant-level z across pathways; white flipped; dorsal=hit_rate raw] ===",
  capture.output(summary(combined_model_signed_corridor))
)

writeLines(summary_text, con = output_txt)
cat("Model results saved:", output_txt, "\n")
