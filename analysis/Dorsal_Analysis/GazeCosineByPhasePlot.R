library(dplyr)
library(ggplot2)

# Argument 1: gaze cosine CSV. Argument 2: output path.
default_input <- "data/processed/gaze_target_cosine_occlusion.csv"
default_output <- "outputs/Gaze_cos/gaze_cosine_by_phase.png"

metric_col <- "cos_mean"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output

df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = tolower(phase),
    participant = as.character(participant)
  ) %>%
  filter(phase %in% c("perception", "explicit", "implicit")) %>%
  filter(is.finite(.data[[metric_col]]))

df_sum <- df %>%
  group_by(participant, phase) %>%
  summarize(value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(value))

keep_ids <- df_sum %>%
  group_by(participant) %>%
  summarize(n_phase = n_distinct(phase), .groups = "drop") %>%
  filter(n_phase == 3) %>%
  pull(participant)

df_sum <- df_sum %>%
  filter(participant %in% keep_ids) %>%
  mutate(
    phase = factor(phase, levels = c("perception", "explicit", "implicit"))
  )

p <- ggplot(df_sum, aes(x = phase, y = value, group = participant)) +
  geom_line(alpha = 0.3, color = "#7A7A7A") +
  geom_point(size = 1.8, alpha = 0.8, color = "#2C7FB8") +
  labs(
    title = "Dorsal cosine by phase",
    x = "Phase",
    y = "Cosine (mean)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(output_path, p, width = 6.5, height = 4.8, dpi = 300)
cat("Saved:", output_path, "\n")
