library(dplyr)

# ==== Input/Output ====
# Argument 1: RDS directory or a single RDS file. Argument 2: optional output directory.
default_input <- "data/processed/RDS_Data"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  input_path <- normalizePath(args[1], mustWork = TRUE)
} else {
  input_path <- default_input
}
if (length(args) >= 2) {
  output_dir <- normalizePath(args[2], mustWork = FALSE)
} else {
  output_dir <- default_output_dir
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

trial_output <- file.path(output_dir, "ventral_trial_responses.csv")

if (dir.exists(input_path)) {
  rds_files <- list.files(input_path, pattern = "^Data_.*\\.rds$", full.names = TRUE, ignore.case = TRUE)
  rds_files <- sort(rds_files)
} else if (file.exists(input_path)) {
  rds_files <- input_path
} else {
  stop("RDS file or directory not found: ", input_path)
}

if (!length(rds_files)) {
  stop("No RDS files found in: ", input_path)
}

# Extract the stimulus color label from each video filename.
extract_stimulus <- function(video_name) {
  if (is.na(video_name)) return(NA_character_)
  if (grepl("white", video_name, ignore.case = TRUE)) return("white")
  if (grepl("black", video_name, ignore.case = TRUE)) return("black")
  NA_character_
}

# Extract the trial suffix embedded at the end of each video filename.
extract_suffix <- function(video_name) {
  if (is.na(video_name)) return(NA_character_)
  base <- tools::file_path_sans_ext(basename(video_name))
  suffix <- sub(".*-(\\d+)$", "\\1", base)
  if (identical(suffix, base)) return(NA_character_)
  suffix
}

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

trial_results <- list()

for (rds_file in rds_files) {
  cat("Processing RDS:", rds_file, "\n")
  data <- readRDS(rds_file)
  if (!nrow(data)) {
    cat("Warning: empty RDS, skipping:", rds_file, "\n")
    next
  }

  required_cols <- c("TRIAL_INDEX", "Time", "Pupil", "phase", "video_file_p", "video_file_i", "video_file_e")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols)) {
    cat("Warning: missing required columns, skipping:", rds_file, "-", paste(missing_cols, collapse = ", "), "\n")
    next
  }

  participant_name <- basename(rds_file)
  participant_name <- sub("^Data_", "", participant_name)
  participant_name <- sub("\\.[Rr][Dd][Ss]$", "", participant_name)
  participant_name <- pad_id(participant_name)

  data <- data %>%
    mutate(
      video_file_mapped = case_when(
        phase == "phase_p" ~ video_file_p,
        phase == "phase_i" ~ video_file_i,
        phase == "phase_e" ~ video_file_e,
        TRUE ~ NA_character_
      ),
      stimulus = vapply(video_file_mapped, extract_stimulus, character(1)),
      suffix = vapply(video_file_mapped, extract_suffix, character(1))
    ) %>%
    filter(phase %in% c("phase_p", "phase_i", "phase_e")) %>%
    filter(!is.na(video_file_mapped) & video_file_mapped != ".") %>%
    filter(!is.na(stimulus)) %>%
    filter(phase == "phase_p" | suffix == "1")

  trial_base <- data %>%
    filter(Time >= 0, Time <= 7000) %>%
    group_by(participant = participant_name, phase, stimulus, video_file_mapped, TRIAL_INDEX) %>%
    summarize(baseline = mean(Pupil, na.rm = TRUE), .groups = "drop")

  trial_roi <- data %>%
    filter(Time >= 8000, Time <= 12000) %>%
    group_by(participant = participant_name, phase, stimulus, video_file_mapped, TRIAL_INDEX) %>%
    summarize(roi = mean(Pupil, na.rm = TRUE), .groups = "drop")

  trial_all <- full_join(
    trial_base, trial_roi,
    by = c("participant", "phase", "stimulus", "video_file_mapped", "TRIAL_INDEX")
  ) %>%
    mutate(pupil_response = roi - baseline)

  trial_results[[length(trial_results) + 1]] <- trial_all
}

if (!length(trial_results)) {
  stop("No usable results were generated.")
}

trial_df <- bind_rows(trial_results)

# ==== Participant-Level Summary ====
stim_means <- trial_df %>%
  group_by(participant, phase, stimulus) %>%
  summarize(pupil_response_mean = mean(pupil_response, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = stimulus,
    values_from = pupil_response_mean,
    names_prefix = "pupil_response_mean_"
  ) %>%
  mutate(pupil_response_diff_bw = pupil_response_mean_black - pupil_response_mean_white)

trial_df <- trial_df %>%
  left_join(stim_means, by = c("participant", "phase"))

write.csv(trial_df, trial_output, row.names = FALSE)
cat("Saved:", trial_output, "\n")
