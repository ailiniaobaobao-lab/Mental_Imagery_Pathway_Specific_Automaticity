library(dplyr)
library(tidyr)
library(lme4)

# ==== Input/Output ====
# Argument 1: ventral trial CSV. Argument 2: dorsal trial CSV. Argument 3: optional output directory.
# Argument 4: top_pct (default 0.30).
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_dorsal <- "data/processed/dorsal_trial_errors.csv"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
dorsal_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_dorsal
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
top_pct <- if (length(args) >= 4) as.numeric(args[4]) else 0.30

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_slopes <- file.path(output_dir, "automaticity_onoff_slopes.csv")
output_summary <- file.path(output_dir, "automaticity_onoff_summary.csv")

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

# Filter and reformat the standardized data for phase-slope modeling.
prepare_model_data <- function(df) {
  df %>%
    filter(phase %in% c("explicit", "implicit")) %>%
    mutate(
      participant = factor(participant),
      phase = factor(phase, levels = c("explicit", "implicit")),
      stimulus = factor(stimulus, levels = c("white", "black"))
    ) %>%
    group_by(participant) %>%
    filter(sum(phase == "explicit") >= 2, sum(phase == "implicit") >= 2) %>%
    ungroup()
}

fit_phase_slope <- function(df) {
  lmer(dv_z ~ phase + stimulus + (1 + phase | participant), data = df, REML = FALSE)
}

# Extract participant-specific random slopes for the target phase contrast.
extract_phase_slopes <- function(model, phase_level = "implicit") {
  re <- ranef(model)$participant
  phase_name <- paste0("phase", phase_level)
  if (!phase_name %in% names(fixef(model))) {
    stop("Missing fixed effect: ", phase_name)
  }
  if (!phase_name %in% names(re)) {
    stop("Missing random slope: ", phase_name)
  }
  slope <- fixef(model)[[phase_name]] + re[[phase_name]]
  tibble(
    participant = rownames(re),
    slope = as.numeric(slope)
  )
}

flag_top_pct <- function(x, top_pct) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  thr <- stats::quantile(abs(x), probs = 1 - top_pct, na.rm = TRUE, names = FALSE, type = 7)
  abs(x) >= thr
}

# ==== Ventral: white sign-flip + perception standardization ====
ventral <- read.csv(ventral_path, stringsAsFactors = FALSE)
ventral_clean <- ventral %>%
  mutate(
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    participant = pad_id(participant),
    dv_raw = ifelse(stimulus == "white", -pupil_response, pupil_response)
  )
ventral_z <- add_perception_z(ventral_clean, "dv_raw")
ventral_model_df <- prepare_model_data(ventral_z)
ventral_model <- fit_phase_slope(ventral_model_df)
ventral_slopes <- extract_phase_slopes(ventral_model) %>%
  rename(slope_ventral = slope)

# ==== Dorsal: strength + perception standardization ====
strength_eps <- 1e-6
dorsal <- read.csv(dorsal_path, stringsAsFactors = FALSE)
dorsal_clean <- dorsal %>%
  mutate(
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    participant = pad_id(participant),
    strength = 1 / (pmax(traj_rmse, 0) + strength_eps)
  )
dorsal_z <- add_perception_z(dorsal_clean, "strength")
dorsal_model_df <- prepare_model_data(dorsal_z)
dorsal_model <- fit_phase_slope(dorsal_model_df)
dorsal_slopes <- extract_phase_slopes(dorsal_model) %>%
  rename(slope_dorsal = slope)

slopes <- full_join(ventral_slopes, dorsal_slopes, by = "participant") %>%
  mutate(
    on_ventral = flag_top_pct(slope_ventral, top_pct),
    on_dorsal = flag_top_pct(slope_dorsal, top_pct)
  )

valid <- slopes %>% filter(!is.na(on_ventral), !is.na(on_dorsal))
a <- sum(valid$on_ventral & valid$on_dorsal)
b <- sum(valid$on_ventral & !valid$on_dorsal)
c <- sum(!valid$on_ventral & valid$on_dorsal)
d <- sum(!valid$on_ventral & !valid$on_dorsal)
den <- (a + b) * (c + d) * (a + c) * (b + d)
phi <- if (den == 0) NA_real_ else (a * d - b * c) / sqrt(den)

summary_df <- tibble(
  top_pct = top_pct,
  n_participants = nrow(valid),
  cor_slope = cor(valid$slope_dorsal, valid$slope_ventral, use = "complete.obs"),
  overlap_on = a,
  overlap_total = nrow(valid),
  phi = phi
)

write.csv(slopes, output_slopes, row.names = FALSE)
write.csv(summary_df, output_summary, row.names = FALSE)

cat("Saved:", output_slopes, "\n")
cat("Saved:", output_summary, "\n")
