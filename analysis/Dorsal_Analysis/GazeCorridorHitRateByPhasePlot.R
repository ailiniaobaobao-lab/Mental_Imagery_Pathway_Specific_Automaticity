library(dplyr)
library(ggplot2)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}

# Argument 1: corridor trial CSV. Argument 2: output path.
default_input <- "data/processed/gaze_corridor_hit_rate_trials.csv"
default_output <- "outputs/Gaze_corrodor/gaze_corridor_hit_rate_diff_ie.png"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output

# Load corridor hit rates and compute participant-level means by phase.
df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = tolower(phase),
    participant = as.character(participant)
  ) %>%
  filter(phase %in% c("perception", "explicit", "implicit")) %>%
  filter(is.finite(hit_rate))

df_sum <- df %>%
  group_by(participant, phase) %>%
  summarize(hit_rate = mean(hit_rate, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(hit_rate))

wide <- df_sum %>%
  filter(phase %in% c("explicit", "implicit")) %>%
  tidyr::pivot_wider(names_from = phase, values_from = hit_rate) %>%
  filter(is.finite(explicit), is.finite(implicit)) %>%
  mutate(diff_ie = implicit - explicit)

# Visualize the implicit-minus-explicit difference distribution across participants.
p <- ggplot(wide, aes(x = "", y = diff_ie)) +
  geom_hline(yintercept = 0, color = "#8E8E8E", linewidth = 0.6) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = "#9EC9E2", color = "#2C7FB8") +
  geom_jitter(width = 0.08, height = 0, color = "#2C7FB8", size = 1.6, alpha = 0.7) +
  labs(
    title = "Corridor hit rate: implicit - explicit",
    x = "",
    y = "Hit rate difference"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(output_path, p, width = 6.5, height = 4.8, dpi = 300)
cat("Saved:", output_path, "\n")
