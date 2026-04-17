library(dplyr)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}
library(emmeans)
library(ggplot2)
library(readxl)

# Argument 1: ventral trial CSV. Argument 2: participants xlsx. Argument 3: output directory.
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_info <- "data/metadata/Participants_info.xlsx"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_txt <- file.path(output_dir, "ventral_group_dorsal_mixed_model.txt")
output_emm <- file.path(output_dir, "ventral_group_dorsal_emmeans.csv")
output_plot <- file.path(output_dir, "ventral_group_dorsal_emmeans.png")

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

# Load dorsal subgroup labels so ventral responses can be modeled by dorsal grouping.
info <- readxl::read_excel(info_path)
col_h <- names(info)[8]
group_map <- info %>%
  transmute(
    participant = pad_id(.data[["subject#"]]),
    group_dorsal = as.character(.data[[col_h]])
  ) %>%
  filter(group_dorsal %in% c("B", "N", "E"))

# Load ventral trial data and flip white trials so the sign reflects the same direction.
ventral_raw <- read.csv(ventral_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    trial_id = sub(".*-(\\d+)\\..*$", "\\1", as.character(video_file_mapped)),
    dv_raw = ifelse(stimulus == "white", -pupil_response, pupil_response)
  )

# Standardize ventral responses relative to perception and collapse to trial-level summaries.
ventral_trial <- add_perception_z(ventral_raw, "dv_raw") %>%
  filter(!is.na(trial_id)) %>%
  group_by(participant, phase, stimulus, trial_id) %>%
  summarize(ventral_value = mean(dv_z, na.rm = TRUE), .groups = "drop") %>%
  inner_join(group_map, by = "participant") %>%
  mutate(
    participant = factor(participant),
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    stimulus = factor(stimulus, levels = c("white", "black")),
    group_dorsal = factor(group_dorsal, levels = c("B", "N", "E"))
  )

# Fit the mixed model, then export emmeans and a matching summary plot.
model <- lmer(ventral_value ~ phase * group_dorsal + stimulus + (1 | participant),
              data = ventral_trial, REML = FALSE)

summary_text <- c(
  "=== Ventral mixed model with dorsal group: ventral_value ~ phase * group_dorsal + stimulus + (1|participant) ===",
  capture.output(summary(model))
)

writeLines(summary_text, con = output_txt)

emm_df <- as.data.frame(emmeans(model, ~ phase * group_dorsal)) %>%
  mutate(section = "emmeans_phase_group")
write.csv(emm_df, output_emm, row.names = FALSE)

plot_df <- emm_df %>%
  mutate(
    phase = factor(phase, levels = c("perception", "explicit", "implicit")),
    group_dorsal = factor(group_dorsal, levels = c("B", "E", "N"))
  )

p <- ggplot(plot_df, aes(x = phase, y = emmean, color = group_dorsal, group = group_dorsal)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.08, linewidth = 0.6) +
  labs(
    title = "Ventral value by phase (dorsal group)",
    x = "Phase",
    y = "Ventral value (perception-z)",
    color = "Group"
  ) +
  theme_minimal(base_size = 12)

ggsave(output_plot, p, width = 6.5, height = 4.5, dpi = 300)

cat("Saved:", output_txt, "\n")
cat("Saved:", output_emm, "\n")
cat("Saved:", output_plot, "\n")
