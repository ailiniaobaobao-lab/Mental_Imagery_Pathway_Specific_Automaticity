library(dplyr)
library(tidyr)

# ==== Input/Output ====
# Argument 1: ventral trial CSV. Argument 2: dorsal trial CSV. Argument 3: optional output directory.
# Argument 4: number of split iterations (default 200). Argument 5: minimum samples per phase-half (default 1).
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_dorsal <- "data/processed/dorsal_trial_errors.csv"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
dorsal_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_dorsal
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
split_n <- if (length(args) >= 4) as.integer(args[4]) else 200L
min_per_phase_half <- if (length(args) >= 5) as.integer(args[5]) else 1L

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_runs <- file.path(output_dir, "automaticity_split_half_runs.csv")
output_summary <- file.path(output_dir, "automaticity_split_half_summary.csv")

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

# Randomly split each participant-phase dataset into two balanced halves.
split_half_assign <- function(df) {
  df %>%
    group_by(participant, phase) %>%
    mutate(
      split_id = {
        n <- n()
        idx <- sample.int(n)
        cut <- ceiling(n / 2)
        rep(c(1, 2), times = c(cut, n - cut))[idx]
      }
    ) %>%
    ungroup()
}

estimate_slope <- function(df) {
  if (!nrow(df)) return(NA_real_)
  phases <- unique(df$phase)
  if (!all(c("explicit", "implicit") %in% phases)) return(NA_real_)
  if (n_distinct(df$stimulus) >= 2) {
    fit <- try(stats::lm(dv_z ~ phase + stimulus, data = df), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      coef_name <- "phaseimplicit"
      if (coef_name %in% names(coef(fit))) return(as.numeric(coef(fit)[[coef_name]]))
    }
  }
  mean(df$dv_z[df$phase == "implicit"], na.rm = TRUE) -
    mean(df$dv_z[df$phase == "explicit"], na.rm = TRUE)
}

# Estimate split-half reliability across repeated random half assignments.
compute_split_reliability <- function(df, split_n, min_per_phase_half, pathway_label) {
  results <- vector("list", split_n)
  for (i in seq_len(split_n)) {
    split_df <- split_half_assign(df)
    slopes <- split_df %>%
      group_by(participant, split_id, phase) %>%
      filter(n() >= min_per_phase_half) %>%
      ungroup() %>%
      group_by(participant, split_id) %>%
      group_modify(~tibble(slope = estimate_slope(.x))) %>%
      ungroup()

    wide <- slopes %>%
      pivot_wider(names_from = split_id, values_from = slope, names_prefix = "half_")

    r <- cor(wide$half_1, wide$half_2, use = "complete.obs")
    results[[i]] <- tibble(
      iteration = i,
      pathway = pathway_label,
      r = r,
      n_participants = sum(complete.cases(wide$half_1, wide$half_2))
    )
  }
  bind_rows(results)
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
ventral_z <- add_perception_z(ventral_clean, "dv_raw") %>%
  filter(phase %in% c("explicit", "implicit"))

ventral_runs <- compute_split_reliability(ventral_z, split_n, min_per_phase_half, "ventral")

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
dorsal_z <- add_perception_z(dorsal_clean, "strength") %>%
  filter(phase %in% c("explicit", "implicit"))

dorsal_runs <- compute_split_reliability(dorsal_z, split_n, min_per_phase_half, "dorsal")

all_runs <- bind_rows(ventral_runs, dorsal_runs)

summary_df <- all_runs %>%
  group_by(pathway) %>%
  summarize(
    split_n = split_n,
    r_mean = mean(r, na.rm = TRUE),
    r_median = median(r, na.rm = TRUE),
    ci_low = quantile(r, 0.025, na.rm = TRUE, names = FALSE, type = 7),
    ci_high = quantile(r, 0.975, na.rm = TRUE, names = FALSE, type = 7),
    n_participants = round(mean(n_participants, na.rm = TRUE)),
    .groups = "drop"
  )

write.csv(all_runs, output_runs, row.names = FALSE)
write.csv(summary_df, output_summary, row.names = FALSE)

cat("Saved:", output_runs, "\n")
cat("Saved:", output_summary, "\n")
