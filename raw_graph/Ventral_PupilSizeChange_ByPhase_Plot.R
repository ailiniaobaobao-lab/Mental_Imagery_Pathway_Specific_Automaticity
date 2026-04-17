library(dplyr)
library(tidyr)
library(ggplot2)

# ==== Input/Output ====
# Argument 1: ventral trial CSV. Argument 2: optional output directory.
default_input <- "data/processed/ventral_trial_responses.csv"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_dir <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

data <- read.csv(input_path, stringsAsFactors = FALSE)
required_cols <- c("participant", "phase", "stimulus", "pupil_response")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols)) {
  stop("Input file is missing columns: ", paste(missing_cols, collapse = ", "))
}

exclude_participants <- c(11191401, 11191503, 11210905, 11241311, 11241210)

data <- data %>%
  mutate(
    phase = case_when(
      phase == "phase_i" ~ "implicit",
      phase == "phase_e" ~ "explicit",
      TRUE ~ phase
    ),
    phase = factor(phase, levels = c("explicit", "implicit")),
    stimulus = factor(stimulus, levels = c("white", "black"))
  )

delta_df <- data %>%
  group_by(participant, phase, stimulus) %>%
  summarize(mean_response = mean(pupil_response, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = stimulus,
    values_from = mean_response,
    id_cols = c(participant, phase)
  ) %>%
  mutate(delta = black - white) %>%
  filter(!is.na(delta), !participant %in% exclude_participants)

participant_order <- delta_df %>%
  distinct(participant) %>%
  arrange(participant) %>%
  pull(participant)

delta_df <- delta_df %>%
  mutate(participant = factor(participant, levels = participant_order))

p <- ggplot(delta_df, aes(x = participant, y = delta, color = phase)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
  labs(
    title = "Ventral Δ (black − white): Explicit vs Implicit",
    x = "Participant",
    y = "Δ pupil response (mm)",
    color = "Phase"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

output_file <- file.path(output_dir, "ventral_delta_combined.png")
ggsave(output_file, plot = p, width = 12, height = 5.5, dpi = 300)
cat("Saved:", output_file, "\n")
