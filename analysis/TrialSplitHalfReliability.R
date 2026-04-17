library(dplyr)
library(tidyr)

# Argument 1: ventral trial CSV. Argument 2: corridor trial CSV. Argument 3: output directory.
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_corridor <- "data/processed/gaze_corridor_hit_rate_trials.csv"
default_output_dir <- "outputs/Trial_decoupling"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
corridor_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_corridor
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_txt <- file.path(output_dir, "trial_split_half_reliability.txt")
output_csv <- file.path(output_dir, "trial_split_half_reliability_samples.csv")

n_iter <- 200
phases_keep <- c("perception", "explicit", "implicit")

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
    TRUE ~ tolower(as.character(x))
  )
}

parse_suffix <- function(x) {
  out <- sub(".*-(\\d+)\\..*$", "\\1", as.character(x))
  out <- ifelse(out == as.character(x), NA_character_, out)
  out
}

assign_halves_stratified <- function(df, stratum_col) {
  if (!nrow(df)) return(df)
  df <- df %>%
    group_by(.data[[stratum_col]]) %>%
    mutate(
      rand = runif(n()),
      order = rank(rand, ties.method = "random"),
      half = ifelse(order <= floor(n() / 2), "A", "B")
    ) %>%
    ungroup() %>%
    select(-rand, -order)
  df
}

split_half_corr <- function(df, metric_col, n_iter, phases, stratum_col) {
  results <- list()
  for (phase_label in phases) {
    phase_df <- df %>% filter(phase == phase_label)
    if (!nrow(phase_df)) next
    participants <- sort(unique(phase_df$participant))
    for (iter in seq_len(n_iter)) {
      half_rows <- list()
      for (pid in participants) {
        sub <- phase_df %>% filter(participant == pid)
        sub <- assign_halves_stratified(sub, stratum_col)
        half1 <- sub %>% filter(half == "A") %>% summarize(v = mean(.data[[metric_col]], na.rm = TRUE))
        half2 <- sub %>% filter(half == "B") %>% summarize(v = mean(.data[[metric_col]], na.rm = TRUE))
        if (!is.finite(half1$v) || !is.finite(half2$v)) next
        half_rows[[length(half_rows) + 1]] <- data.frame(
          participant = pid,
          half1 = half1$v,
          half2 = half2$v
        )
      }
      if (!length(half_rows)) next
      half_df <- bind_rows(half_rows)
      if (nrow(half_df) < 3) {
        r_val <- NA_real_
        n_val <- nrow(half_df)
      } else {
        r_val <- suppressWarnings(stats::cor(half_df$half1, half_df$half2, use = "complete.obs"))
        n_val <- nrow(half_df)
      }
      results[[length(results) + 1]] <- data.frame(
        phase = phase_label,
        iter = iter,
        r = r_val,
        n = n_val
      )
    }
  }
  bind_rows(results)
}

summarize_reliability <- function(res_df) {
  res_df %>%
    group_by(phase) %>%
    summarize(
      r_mean = mean(r, na.rm = TRUE),
      r_median = median(r, na.rm = TRUE),
      r_ci_low = quantile(r, 0.025, na.rm = TRUE),
      r_ci_high = quantile(r, 0.975, na.rm = TRUE),
      n_mean = mean(n, na.rm = TRUE),
      .groups = "drop"
    )
}

# ==== Ventral: signed pupil_response (perception-z using full baseline) ====
ventral_raw <- read.csv(ventral_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    participant = pad_id(participant),
    trial_id = parse_suffix(video_file_mapped),
    dv_raw = ifelse(stimulus == "white", -pupil_response, pupil_response)
  ) %>%
  filter(phase %in% phases_keep, !is.na(trial_id), !is.na(dv_raw))

ventral_base <- ventral_raw %>%
  filter(phase == "perception") %>%
  group_by(participant, stimulus) %>%
  summarize(
    base_mean = mean(dv_raw, na.rm = TRUE),
    base_sd = sd(dv_raw, na.rm = TRUE),
    .groups = "drop"
  )

ventral <- ventral_raw %>%
  left_join(ventral_base, by = c("participant", "stimulus")) %>%
  mutate(metric = (dv_raw - base_mean) / base_sd) %>%
  filter(is.finite(metric)) %>%
  mutate(
    trial_key = paste(stimulus, trial_id, sep = "_"),
    stratum = trial_key
  )

# ==== Corridor: hit rate (raw) ====
corridor <- read.csv(corridor_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    participant = pad_id(participant),
    trial_id = as.character(trial_id),
    trial_key = paste(stimulus, trial_id, sep = "_"),
    stage = as.character(stage),
    metric = as.numeric(hit_rate),
    stratum = paste(stimulus, stage, trial_id, sep = "_")
  ) %>%
  filter(phase %in% phases_keep, !is.na(trial_id), !is.na(metric))

ventral_res <- split_half_corr(ventral, "metric", n_iter, phases_keep, "stratum") %>%
  mutate(pathway = "ventral")
corridor_res <- split_half_corr(corridor, "metric", n_iter, phases_keep, "stratum") %>%
  mutate(pathway = "dorsal")

all_res <- bind_rows(ventral_res, corridor_res)
write.csv(all_res, output_csv, row.names = FALSE)

ventral_sum <- summarize_reliability(ventral_res) %>% mutate(pathway = "ventral")
corridor_sum <- summarize_reliability(corridor_res) %>% mutate(pathway = "dorsal")
summary_df <- bind_rows(ventral_sum, corridor_sum)

lines <- c(
  "=== Split-half reliability (trial-level, stratified) ===",
  paste0("Iterations: ", n_iter),
  "",
  "Format: r_mean (median) with 95% CI, mean n per iteration"
)

for (pathway_label in c("ventral", "dorsal")) {
  lines <- c(lines, "", paste0("== ", pathway_label, " =="))
  sub <- summary_df %>% filter(pathway == pathway_label)
  for (phase_label in phases_keep) {
    row <- sub %>% filter(phase == phase_label)
    if (!nrow(row)) next
    lines <- c(
      lines,
      paste0(
        phase_label, ": r_mean=", sprintf("%.3f", row$r_mean),
        " (median=", sprintf("%.3f", row$r_median),
        "), 95% CI [", sprintf("%.3f", row$r_ci_low), ", ", sprintf("%.3f", row$r_ci_high),
        "], n_mean=", sprintf("%.1f", row$n_mean)
      )
    )
  }
}

writeLines(lines, con = output_txt)
cat("Saved:", output_txt, "\n")
cat("Saved:", output_csv, "\n")
