library(readxl)
library(dplyr)
library(ggplot2)

# Argument 1: corridor trial CSV. Argument 2: participants xlsx. Argument 3: output path.
default_input <- "data/processed/gaze_corridor_hit_rate_trials.csv"
default_info <- "data/metadata/Participants_info.xlsx"
default_output <- "outputs/Gaze_corrodor/gaze_corridor_hit_rate_by_phase_group_mean.png"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_path <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output

normalize_phase <- function(x) {
  dplyr::case_when(
    x == "phase_p" ~ "perception",
    x == "phase_e" ~ "explicit",
    x == "phase_i" ~ "implicit",
    TRUE ~ as.character(x)
  )
}

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

# Load participant group labels so the group mean plot uses the same subgroup mapping.
info <- readxl::read_excel(info_path)
col_h <- names(info)[8]
group_map <- info %>%
  transmute(
    participant = pad_id(.data[["subject#"]]),
    group = as.character(.data[[col_h]])
  ) %>%
  filter(group %in% c("B", "E", "N"))

# Load corridor hit rates and summarize them within participant and phase.
df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = normalize_phase(phase)
  ) %>%
  filter(is.finite(hit_rate), phase %in% c("perception", "explicit", "implicit"))

df_sum <- df %>%
  group_by(participant, phase) %>%
  summarize(hit_rate = mean(hit_rate, na.rm = TRUE), .groups = "drop") %>%
  inner_join(group_map, by = "participant")

keep_ids <- df_sum %>%
  group_by(participant) %>%
  summarize(n_phase = n_distinct(phase), .groups = "drop") %>%
  filter(n_phase == 3) %>%
  pull(participant)

df_sum <- df_sum %>%
  filter(participant %in% keep_ids) %>%
  mutate(phase = factor(phase, levels = c("perception", "explicit", "implicit")))

# Aggregate participant summaries to group-level means and confidence intervals.
group_stats <- df_sum %>%
  group_by(group, phase) %>%
  summarize(
    mean_hit = mean(hit_rate, na.rm = TRUE),
    sd_hit = sd(hit_rate, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se = ifelse(n > 1, sd_hit / sqrt(n), NA_real_),
    t_crit = ifelse(n > 1, stats::qt(0.975, df = n - 1), NA_real_),
    ci = se * t_crit,
    ymin = mean_hit - ci,
    ymax = mean_hit + ci
  )

phase_labels <- c(
  perception = "Perception",
  explicit = "Instructed Imagery",
  implicit = "Passive Viewing"
)

group_labels <- c(
  B = "B",
  E = "I",
  N = "N"
)

# Plot the subgroup means with 95% confidence intervals across phases.
p <- ggplot(group_stats, aes(x = phase, y = mean_hit, color = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.08, linewidth = 0.6) +
  labs(
    x = "Phase",
    y = "Corridor hit rate (mean ± 95% CI)",
    color = "Group"
  ) +
  scale_x_discrete(labels = phase_labels) +
  scale_color_discrete(labels = group_labels) +
  theme_minimal(base_size = 13) +
  theme(text = element_text(family = "Times New Roman"))

ggsave(output_path, p, width = 6.8, height = 4.6, dpi = 300)
cat("Saved:", output_path, "\n")
