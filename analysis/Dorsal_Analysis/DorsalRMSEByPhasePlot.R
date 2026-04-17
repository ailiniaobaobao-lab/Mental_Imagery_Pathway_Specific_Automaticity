library(dplyr)
library(ggplot2)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}

# Argument 1: dorsal trial CSV. Argument 2: output path.
default_input <- "data/processed/dorsal_trial_errors.csv"
default_output <- "outputs/PupilSize_Raw/dorsal_rmse_by_phase.png"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output

# Load dorsal trial data and summarize RMSE within participant and phase.
df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = tolower(phase),
    participant = as.character(participant),
    traj_rmse = as.numeric(traj_rmse)
  ) %>%
  filter(phase %in% c("perception", "explicit", "implicit")) %>%
  filter(is.finite(traj_rmse))

df_sum <- df %>%
  group_by(participant, phase) %>%
  summarize(traj_rmse = mean(traj_rmse, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(traj_rmse))

keep_ids <- df_sum %>%
  group_by(participant) %>%
  summarize(n_phase = n_distinct(phase), .groups = "drop") %>%
  filter(n_phase == 3) %>%
  pull(participant)

df_sum <- df_sum %>%
  filter(participant %in% keep_ids) %>%
  mutate(
    phase = factor(phase, levels = c("perception", "explicit", "implicit")),
    phase_num = as.numeric(phase)
  )

# Test the linear phase trend and annotate the resulting participant trajectory plot.
model <- lmer(traj_rmse ~ phase_num + (1 | participant), data = df_sum, REML = FALSE)
coef_df <- summary(model)$coefficients
if ("phase_num" %in% rownames(coef_df) && "Pr(>|t|)" %in% colnames(coef_df)) {
  p_val <- coef_df["phase_num", "Pr(>|t|)"]
} else {
  p_val <- NA_real_
}

sig_label <- if (is.na(p_val)) {
  "p = NA"
} else if (p_val < 0.001) {
  sprintf("trend p = %.3g ***", p_val)
} else if (p_val < 0.01) {
  sprintf("trend p = %.3g **", p_val)
} else if (p_val < 0.05) {
  sprintf("trend p = %.3g *", p_val)
} else {
  sprintf("trend p = %.3g", p_val)
}

p <- ggplot(df_sum, aes(x = phase, y = traj_rmse, group = participant)) +
  geom_line(alpha = 0.3, color = "#7A7A7A") +
  geom_point(size = 1.8, alpha = 0.8, color = "#2C7FB8") +
  annotate(
    "text",
    x = 3,
    y = Inf,
    label = sig_label,
    hjust = 1.02,
    vjust = 1.2,
    size = 3.5,
    color = "grey20"
  ) +
  labs(
    title = "Dorsal RMSE by phase",
    x = "Phase",
    y = "RMSE"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(output_path, p, width = 6.5, height = 4.8, dpi = 300)
cat("Saved:", output_path, "\n")
