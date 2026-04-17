library(dplyr) # data wrangling
library(ggplot2) # plotting

# ==== Input/Output ====
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_dorsal <- "data/processed/dorsal_trial_errors.csv"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
dorsal_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_dorsal
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

out_exp_imp <- file.path(output_dir, "combined_re_strength_explicit_implicit.png")
out_exp_imp_per <- file.path(output_dir, "combined_re_strength_explicit_implicit_perception.png")

strength_eps <- 1e-6
exclude_participants <- c(11191401, 11191503, 11210905, 11241311, 11241210)

# Load ventral and dorsal trial tables used to build the combined strength plots.
ventral <- read.csv(ventral_path, stringsAsFactors = FALSE)
dorsal <- read.csv(dorsal_path, stringsAsFactors = FALSE)

ventral_df <- ventral %>%
  mutate(
    phase = case_when(
      phase == "phase_i" ~ "implicit",
      phase == "phase_e" ~ "explicit",
      phase == "phase_p" ~ "perception",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pupil_response), !is.na(phase), !is.na(participant), !is.na(stimulus)) %>%
  filter(!participant %in% exclude_participants) %>%
  group_by(participant, phase, stimulus) %>%
  summarize(mean_response = mean(pupil_response, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = stimulus, values_from = mean_response) %>%
  mutate(delta_bw = black - white) %>%
  filter(!is.na(delta_bw)) %>%
  transmute(
    participant,
    phase,
    pathway = "ventral",
    dv_raw = delta_bw
  )

dorsal_df <- dorsal %>%
  mutate(
    phase = as.character(phase),
    strength_re = 1 / (pmax(radial_error, 0) + strength_eps)
  ) %>%
  filter(!is.na(strength_re), !is.na(phase), !is.na(participant)) %>%
  filter(!participant %in% exclude_participants) %>%
  group_by(participant, phase) %>%
  summarize(strength_re = mean(strength_re, na.rm = TRUE), .groups = "drop") %>%
  transmute(
    participant,
    phase,
    pathway = "dorsal",
    dv_raw = strength_re
  )

# Standardize both pathways to a shared plotting scale within each pathway.
build_plot_data <- function(phases) {
  combined <- bind_rows(
    ventral_df %>% filter(phase %in% phases),
    dorsal_df %>% filter(phase %in% phases)
  ) %>%
    mutate(
      phase = factor(phase, levels = phases),
      pathway = factor(pathway, levels = c("ventral", "dorsal"))
    ) %>%
    group_by(pathway) %>%
    mutate(dv_z = as.numeric(scale(dv_raw))) %>%
    ungroup()

  participant_means <- combined %>%
    group_by(participant, phase, pathway) %>%
    summarize(dv_z = mean(dv_z, na.rm = TRUE), .groups = "drop")

  summary_df <- participant_means %>%
    group_by(phase, pathway) %>%
    summarize(
      mean_dv = mean(dv_z, na.rm = TRUE),
      se = sd(dv_z, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )

  list(
    participant_means = participant_means,
    summary_df = summary_df,
    phases = phases
  )
}

# Plot participant summaries together with phase-level means for each pathway.
plot_combined <- function(plot_data, title, out_path) {
  participant_means <- plot_data$participant_means
  summary_df <- plot_data$summary_df
  phases <- plot_data$phases

  p <- ggplot() +
    geom_col(
      data = summary_df,
      aes(x = phase, y = mean_dv, fill = phase),
      width = 0.6,
      alpha = 0.6
    ) +
    geom_errorbar(
      data = summary_df,
      aes(x = phase, y = mean_dv, ymin = mean_dv - se, ymax = mean_dv + se),
      width = 0.15,
      color = "#444444"
    ) +
    geom_line(
      data = participant_means,
      aes(x = phase, y = dv_z, group = participant),
      color = "#666666",
      alpha = 0.35
    ) +
    geom_point(
      data = participant_means,
      aes(x = phase, y = dv_z),
      size = 2,
      alpha = 0.8,
      position = position_jitter(width = 0.05, height = 0)
    ) +
    facet_wrap(~pathway, ncol = 2) +
    labs(
      title = title,
      x = NULL,
      y = "Z score (within pathway; ventral=black-white delta, dorsal=RE strength)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(out_path, p, width = 10, height = 5, dpi = 300)
  cat("✅ Saved:", out_path, "\n")
}

# Compare the two pathways first with explicit/implicit only, then with perception included.
data_exp_imp <- build_plot_data(c("explicit", "implicit"))
plot_combined(
  data_exp_imp,
  "Automaticity by pathway (explicit vs implicit)",
  out_exp_imp
)

data_exp_imp_per <- build_plot_data(c("perception", "explicit", "implicit"))
plot_combined(
  data_exp_imp_per,
  "Automaticity by pathway (explicit/implicit/perception)",
  out_exp_imp_per
)
