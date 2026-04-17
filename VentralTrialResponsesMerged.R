library(dplyr)

# Argument 1: input CSV. Argument 2: output HTML.
default_input <- "data/processed/ventral_trial_responses.csv"
default_output <- "outputs/PupilSize_Raw/ventral_trial_responses_merged.html"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

# Load the trial-level ventral table and sort rows for merged HTML display.
df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(participant = pad_id(participant)) %>%
  arrange(participant, phase, TRIAL_INDEX, stimulus)

cols <- c(
  "participant",
  "phase",
  "stimulus",
  "video_file_mapped",
  "TRIAL_INDEX",
  "baseline",
  "roi",
  "pupil_response",
  "pupil_response_mean_black",
  "pupil_response_mean_white",
  "pupil_response_diff_bw"
)
df <- df[, cols]

df$group_id <- paste(df$participant, df$phase, sep = "||")
group_sizes <- table(df$group_id)
seen <- setNames(rep(0L, length(group_sizes)), names(group_sizes))

fmt <- function(x) {
  if (is.numeric(x)) {
    ifelse(is.na(x), "", format(x, digits = 6, nsmall = 2, trim = TRUE))
  } else {
    ifelse(is.na(x), "", as.character(x))
  }
}

escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

merged_cols <- c(
  "participant",
  "phase",
  "pupil_response_mean_black",
  "pupil_response_mean_white",
  "pupil_response_diff_bw"
)

# Merge repeated participant-phase fields into grouped HTML rows for readability.
header <- paste0("<th>", escape_html(cols), "</th>")
rows <- vector("character", nrow(df))

for (i in seq_len(nrow(df))) {
  gid <- df$group_id[i]
  seen[gid] <- seen[gid] + 1L
  is_first <- seen[gid] == 1L
  span <- as.integer(group_sizes[[gid]])
  cells <- character()
  for (col in cols) {
    if (col %in% merged_cols) {
      if (is_first) {
        val <- fmt(df[[col]][i])
        style <- ""
        if (col == "pupil_response_diff_bw") {
          diff_val <- df[[col]][i]
          if (!is.na(diff_val) && diff_val > 0) {
            style <- " style=\"background-color:#c6efce\""
          }
        }
        cells <- c(cells, paste0("<td rowspan=\"", span, "\"", style, ">", escape_html(val), "</td>"))
      }
    } else {
      val <- fmt(df[[col]][i])
      cells <- c(cells, paste0("<td>", escape_html(val), "</td>"))
    }
  }
  rows[i] <- paste0("<tr>", paste(cells, collapse = ""), "</tr>")
}

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset=\"utf-8\">",
  "<style>table{border-collapse:collapse}th,td{border:1px solid #ccc;padding:6px 8px;font-family:Arial, sans-serif;font-size:12px;vertical-align:top}</style>",
  "</head><body>",
  "<table><thead><tr>", paste(header, collapse = ""), "</tr></thead>",
  "<tbody>", paste(rows, collapse = ""), "</tbody></table>",
  "</body></html>"
)

# Write the merged table as a standalone HTML report.
writeLines(html, output_path)
cat("Saved:", output_path, "\n")
