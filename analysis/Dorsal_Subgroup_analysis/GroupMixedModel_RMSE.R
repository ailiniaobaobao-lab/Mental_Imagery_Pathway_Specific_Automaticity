library(dplyr)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}
library(readxl)

# Argument 1: dorsal trial CSV. Argument 2: participants xlsx. Argument 3: output directory.
default_dorsal <- "data/processed/dorsal_trial_errors.csv"
default_info <- "data/metadata/Participants_info.xlsx"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
dorsal_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_dorsal
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_txt <- file.path(output_dir, "dorsal_group_mixed_model.txt")

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

normalize_phase <- function(x) {
  dplyr::case_when(
    x == "phase_i" ~ "implicit",
    x == "phase_e" ~ "explicit",
    x == "phase_p" ~ "perception",
    TRUE ~ as.character(x)
  )
}

# Load participant group labels so the dorsal model can test group-by-phase effects.
info <- readxl::read_excel(info_path)
col_h <- names(info)[8]
group_map <- info %>%
  transmute(
    participant = pad_id(.data[["subject#"]]),
    group = as.character(.data[[col_h]])
  ) %>%
  filter(group %in% c("B", "N", "E"))

# Load dorsal trial data and convert RMSE into a strength-style metric.
strength_eps <- 1e-6
dorsal_raw <- read.csv(dorsal_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = normalize_phase(phase),
    stimulus = factor(stimulus, levels = c("white", "black")),
    strength_rmse = 1 / (pmax(traj_rmse, 0) + strength_eps)
  ) %>%
  filter(!is.na(strength_rmse), !is.na(participant), !is.na(phase), !is.na(stimulus)) %>%
  inner_join(group_map, by = "participant")

# Use perception trials to standardize each participant and stimulus combination.
perc_stats <- dorsal_raw %>%
  filter(phase == "perception") %>%
  group_by(participant, stimulus) %>%
  summarize(
    base_mean = mean(strength_rmse, na.rm = TRUE),
    base_sd = sd(strength_rmse, na.rm = TRUE),
    .groups = "drop"
  )

# Build the final analysis table and fit the group-by-phase mixed model.
dorsal <- dorsal_raw %>%
  left_join(perc_stats, by = c("participant", "stimulus")) %>%
  mutate(dv_z = (strength_rmse - base_mean) / base_sd) %>%
  filter(is.finite(dv_z)) %>%
  mutate(
    participant = factor(participant),
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    group = factor(group, levels = c("B", "N", "E"))
  )

model <- lmer(dv_z ~ phase * group + stimulus + (1 | participant), data = dorsal, REML = FALSE)

# Export the model summary for downstream reporting.
summary_text <- c(
  "=== Dorsal group mixed model: dv_z (perception-z) ~ phase * group + stimulus + (1|participant) ===",
  capture.output(summary(model))
)

writeLines(summary_text, con = output_txt)
cat("Saved:", output_txt, "\n")
