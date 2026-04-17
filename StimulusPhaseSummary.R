library(readxl)
library(dplyr)
library(tidyr)

# ==== Input/Output ====
# Argument 1: ventral trial CSV. Argument 2: dorsal trial CSV. Argument 3: participant info xlsx. Argument 4: output directory.
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_dorsal <- "data/processed/dorsal_trial_errors.csv"
default_participants <- "data/metadata/Participants_info.xlsx"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
dorsal_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_dorsal
participants_path <- if (length(args) >= 3) normalizePath(args[3], mustWork = TRUE) else default_participants
output_dir <- if (length(args) >= 4) normalizePath(args[4], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

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

# Export a styled HTML table for the phase-by-stimulus summary output.
write_html_table <- function(df, out_path, diff_col) {
  esc <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  fmt <- function(x) {
    if (is.numeric(x)) {
      ifelse(is.na(x), "", format(x, digits = 6, nsmall = 2, trim = TRUE))
    } else {
      ifelse(is.na(x), "", as.character(x))
    }
  }

  headers <- paste0("<th>", esc(names(df)), "</th>")
  rows <- apply(df, 1, function(row) {
    cells <- mapply(function(val, col) {
      val_fmt <- fmt(val)
      style <- ""
      if (col == diff_col && !is.na(as.numeric(val)) && as.numeric(val) > 0) {
        style <- " style=\"background-color:#c6efce\""
      }
      paste0("<td", style, ">", esc(val_fmt), "</td>")
    }, row, names(df), SIMPLIFY = TRUE, USE.NAMES = FALSE)
    paste0("<tr>", paste(cells, collapse = ""), "</tr>")
  })

  html <- paste0(
    "<!DOCTYPE html><html><head><meta charset=\"utf-8\">",
    "<style>table{border-collapse:collapse}th,td{border:1px solid #ccc;padding:6px 8px;font-family:Arial, sans-serif;font-size:12px}</style>",
    "</head><body>",
    "<table><thead><tr>", paste(headers, collapse = ""), "</tr></thead>",
    "<tbody>", paste(rows, collapse = ""), "</tbody></table>",
    "</body></html>"
  )
  writeLines(html, out_path)
}

participants_info <- readxl::read_excel(participants_path)
participant_ids <- participants_info %>%
  transmute(participant = pad_id(`subject#`)) %>%
  filter(!is.na(participant)) %>%
  pull(participant)

# ==== Ventral summary ====
ventral <- read.csv(ventral_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    participant = pad_id(participant)
  ) %>%
  filter(participant %in% participant_ids, phase %in% c("perception", "explicit", "implicit"))

ventral_summary <- ventral %>%
  group_by(participant, phase, stimulus) %>%
  summarize(mean_pupil_response = mean(pupil_response, na.rm = TRUE), n_trials = n(), .groups = "drop") %>%
  pivot_wider(names_from = stimulus, values_from = mean_pupil_response, names_prefix = "mean_") %>%
  mutate(diff_bw = mean_black - mean_white) %>%
  arrange(participant, phase)

ventral_csv <- file.path(output_dir, "ventral_phase_stimulus_mean.csv")
ventral_html <- file.path(output_dir, "ventral_phase_stimulus_mean.html")
write.csv(ventral_summary, ventral_csv, row.names = FALSE)
write_html_table(ventral_summary, ventral_html, "diff_bw")

# ==== Dorsal summary ====
dorsal <- read.csv(dorsal_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    participant = pad_id(participant)
  ) %>%
  filter(participant %in% participant_ids, phase %in% c("perception", "explicit", "implicit"))

dorsal_summary <- dorsal %>%
  group_by(participant, phase, stimulus) %>%
  summarize(mean_traj_rmse = mean(traj_rmse, na.rm = TRUE), n_trials = n(), .groups = "drop") %>%
  pivot_wider(names_from = stimulus, values_from = mean_traj_rmse, names_prefix = "mean_") %>%
  mutate(diff_bw = mean_black - mean_white) %>%
  arrange(participant, phase)

dorsal_csv <- file.path(output_dir, "dorsal_phase_stimulus_mean.csv")
dorsal_html <- file.path(output_dir, "dorsal_phase_stimulus_mean.html")
write.csv(dorsal_summary, dorsal_csv, row.names = FALSE)
write_html_table(dorsal_summary, dorsal_html, "diff_bw")

cat("Saved:", ventral_csv, "\n")
cat("Saved:", ventral_html, "\n")
cat("Saved:", dorsal_csv, "\n")
cat("Saved:", dorsal_html, "\n")
