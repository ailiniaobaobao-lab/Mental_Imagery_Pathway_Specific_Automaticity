library(dplyr)
library(ggplot2)

# Argument 1: dorsal trial CSV. Argument 2: output path.
default_input <- "data/processed/dorsal_trial_errors.csv"
default_output <- "outputs/Gaze_rmse/dorsal_rmse_diff_ie.png"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output

df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = tolower(phase),
    participant = as.character(participant),
    traj_rmse = as.numeric(traj_rmse)
  ) %>%
  filter(phase %in% c("explicit", "implicit")) %>%
  filter(is.finite(traj_rmse))

df_sum <- df %>%
  group_by(participant, phase) %>%
  summarize(traj_rmse = mean(traj_rmse, na.rm = TRUE), .groups = "drop")

wide <- df_sum %>%
  tidyr::pivot_wider(names_from = phase, values_from = traj_rmse) %>%
  filter(is.finite(explicit), is.finite(implicit)) %>%
  mutate(diff_ie = implicit - explicit)

p <- ggplot(wide, aes(x = "", y = diff_ie)) +
  geom_hline(yintercept = 0, color = "#8E8E8E", linewidth = 0.6) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = "#9EC9E2", color = "#2C7FB8") +
  geom_jitter(width = 0.08, height = 0, color = "#2C7FB8", size = 1.6, alpha = 0.7) +
  labs(
    title = "Dorsal RMSE: implicit - explicit",
    x = "",
    y = "RMSE difference"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(output_path, p, width = 4.5, height = 4.5, dpi = 300)
cat("Saved:", output_path, "\n")
