library(dplyr)
library(lme4)
library(ggplot2)
library(emmeans)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}

# Argument 1: ventral trial CSV. Argument 2: output directory.
default_input <- "data/processed/ventral_trial_responses.csv"
default_output_dir <- "outputs/PupilSize"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_dir <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_png <- file.path(output_dir, "ventral_stimulus_effect.png")
output_csv <- file.path(output_dir, "ventral_stimulus_effect_emmeans.csv")

exclude_participants <- c("11191401", "11191503", "11210905", "11241311", "11241210", "03041338")

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

exclude_participants <- pad_id(exclude_participants)

ventral <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = case_when(
      phase == "phase_i" ~ "implicit",
      phase == "phase_e" ~ "explicit",
      phase == "phase_p" ~ "perception",
      TRUE ~ NA_character_
    ),
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    stimulus = factor(stimulus, levels = c("white", "black")),
    participant = pad_id(participant)
  ) %>%
  filter(!is.na(pupil_response), !is.na(participant), !is.na(stimulus), !is.na(phase)) %>%
  filter(!participant %in% exclude_participants)

perc_stats <- ventral %>%
  filter(phase == "perception") %>%
  group_by(participant, stimulus) %>%
  summarize(
    base_mean = mean(pupil_response, na.rm = TRUE),
    base_sd = sd(pupil_response, na.rm = TRUE),
    .groups = "drop"
  )

ventral_z <- ventral %>%
  left_join(perc_stats, by = c("participant", "stimulus")) %>%
  mutate(dv_z = (pupil_response - base_mean) / base_sd) %>%
  filter(is.finite(dv_z)) %>%
  filter(phase %in% c("explicit", "implicit")) %>%
  mutate(participant = factor(participant))

model <- lmer(dv_z ~ phase * stimulus + (1 | participant), data = ventral_z, REML = FALSE)

emm <- emmeans(model, ~ phase * stimulus)
emm_df <- as.data.frame(summary(emm, infer = c(TRUE, TRUE)))
write.csv(emm_df, output_csv, row.names = FALSE)

emm_df <- emm_df %>%
  mutate(
    phase = factor(phase, levels = c("explicit", "implicit")),
    stimulus = factor(stimulus, levels = c("white", "black"))
  )

p <- ggplot(emm_df, aes(x = phase, y = emmean, color = stimulus, group = stimulus)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.08, linewidth = 0.6) +
  scale_color_manual(values = c(white = "#FFD60A", black = "#0B3D91")) +
  labs(
    title = "Ventral stimulus effect (emmeans)",
    x = "Phase",
    y = "Standardized signature strength (perception-z)",
    color = "Stimulus",
    caption = "stimulus beta = -12.040, p = 3.6e-09"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    plot.caption = element_text(size = 9, hjust = 0)
  )

ggsave(output_png, p, width = 6.5, height = 4.5, dpi = 300)
cat("Saved:", output_png, "\n")
cat("Saved:", output_csv, "\n")
