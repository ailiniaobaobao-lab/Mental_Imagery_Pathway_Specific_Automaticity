library(dplyr) # data wrangling
library(tidyr) # data reshaping
library(ggplot2) # plotting
library(lme4) # mixed model

# ==== Input/Output ====
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_dorsal <- "data/processed/dorsal_trial_errors.csv"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
dorsal_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_dorsal
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

out_exp_imp <- file.path(output_dir, "combined_rmse_strength_explicit_implicit.png")
out_exp_imp_per <- file.path(output_dir, "combined_rmse_strength_explicit_implicit_perception.png")
out_exp_imp_delta <- file.path(output_dir, "combined_rmse_strength_explicit_implicit_delta.png")

strength_eps <- 1e-6
exclude_participants <- c(11191401, 11191503, 11210905, 11241311, 11241210)

# Load ventral and dorsal trial tables used to build the cross-pathway comparison plots.
ventral <- read.csv(ventral_path, stringsAsFactors = FALSE)
dorsal <- read.csv(dorsal_path, stringsAsFactors = FALSE)

ventral_df_all <- ventral %>%
  mutate(
    phase = case_when(
      phase == "phase_i" ~ "implicit",
      phase == "phase_e" ~ "explicit",
      phase == "phase_p" ~ "perception",
      TRUE ~ NA_character_
    ),
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    stimulus = factor(stimulus, levels = c("white", "black"))
  ) %>%
  filter(!is.na(pupil_response), !is.na(phase), !is.na(participant), !is.na(stimulus)) %>%
  filter(!participant %in% exclude_participants)

dorsal_df_all <- dorsal %>%
  mutate(
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    stimulus = factor(stimulus, levels = c("white", "black")),
    strength_rmse = 1 / (pmax(traj_rmse, 0) + strength_eps)
  ) %>%
  filter(!is.na(strength_rmse), !is.na(phase), !is.na(participant), !is.na(stimulus)) %>%
  filter(!participant %in% exclude_participants)

# Assemble the shared long-format dataset used across pathway comparison figures.
build_combined_data <- function(phases) {
  ventral_df <- ventral_df_all %>%
    filter(phase %in% phases) %>%
    transmute(
      participant,
      phase,
      pathway = "ventral",
      dv_raw = ifelse(stimulus == "white", -pupil_response, pupil_response)
    )
  dorsal_df <- dorsal_df_all %>%
    filter(phase %in% phases) %>%
    transmute(
      participant,
      phase,
      pathway = "dorsal",
      dv_raw = strength_rmse
    )
  combined <- bind_rows(ventral_df, dorsal_df) %>%
    mutate(
      phase = factor(phase, levels = phases),
      pathway = factor(pathway, levels = c("ventral", "dorsal"))
    ) %>%
    group_by(pathway) %>%
    mutate(dv_z = as.numeric(scale(dv_raw))) %>%
    ungroup()
  combined
}

# Predict fixed-effect means and confidence bands from the mixed model.
predict_fixed <- function(model, phases) {
  newdata <- expand.grid(
    phase = phases,
    pathway = c("ventral", "dorsal"),
    stringsAsFactors = FALSE
  )
  newdata$phase <- factor(newdata$phase, levels = levels(model@frame$phase))
  newdata$pathway <- factor(newdata$pathway, levels = levels(model@frame$pathway))
  terms_obj <- stats::delete.response(terms(model))
  X <- model.matrix(terms_obj, newdata)
  beta <- lme4::fixef(model)
  V <- as.matrix(vcov(model))
  pred <- as.numeric(X %*% beta)
  se <- sqrt(diag(X %*% V %*% t(X)))
  newdata$mean <- pred
  newdata$se <- se
  newdata$lower <- pred - 1.96 * se
  newdata$upper <- pred + 1.96 * se
  newdata
}

# Plot participant trajectories together with the pathway-level model predictions.
plot_combined <- function(combined, title, out_path) {
  phases <- levels(combined$phase)
  model <- lmer(dv_z ~ phase * pathway + (1 | participant), data = combined, REML = FALSE)
  model_means <- predict_fixed(model, phases)

  participant_means <- combined %>%
    group_by(participant, phase, pathway) %>%
    summarize(dv_z = mean(dv_z, na.rm = TRUE), .groups = "drop")

  p <- ggplot() +
    geom_errorbar(
      data = model_means,
      aes(x = phase, y = mean, ymin = lower, ymax = upper),
      width = 0.15,
      color = "#444444"
    ) +
    geom_col(
      data = model_means,
      aes(x = phase, y = mean, fill = phase),
      width = 0.6,
      alpha = 0.6
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
      y = "Z score (within pathway; ventral=flipped pupil response, dorsal=RMSE strength)"
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

# Plot the explicit-versus-implicit pathway comparison first.
data_exp_imp <- build_combined_data(c("explicit", "implicit"))
plot_combined(
  data_exp_imp,
  "Automaticity by pathway (explicit vs implicit)",
  out_exp_imp
)

delta_df <- data_exp_imp %>%
  group_by(participant, pathway, phase) %>%
  summarize(dv_z = mean(dv_z, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = phase, values_from = dv_z) %>%
  mutate(delta = implicit - explicit) %>%
  filter(!is.na(delta))

delta_means <- delta_df %>%
  group_by(pathway) %>%
  summarize(
    mean = mean(delta, na.rm = TRUE),
    se = sd(delta, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_delta <- ggplot() +
  geom_hline(yintercept = 0, color = "#999999") +
  geom_errorbar(
    data = delta_means,
    aes(x = pathway, y = mean, ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
    width = 0.15,
    color = "#444444"
  ) +
  geom_col(
    data = delta_means,
    aes(x = pathway, y = mean, fill = pathway),
    width = 0.6,
    alpha = 0.6
  ) +
  geom_line(
    data = delta_df,
    aes(x = pathway, y = delta, group = participant),
    color = "#666666",
    alpha = 0.35
  ) +
  geom_point(
    data = delta_df,
    aes(x = pathway, y = delta),
    size = 2,
    alpha = 0.8,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  labs(
    title = "Implicit minus explicit (by pathway)",
    x = NULL,
    y = "Delta z-score (implicit - explicit)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(out_exp_imp_delta, p_delta, width = 8, height = 5, dpi = 300)
cat("✅ Saved:", out_exp_imp_delta, "\n")

# Refit the same comparison including perception as the baseline phase.
data_exp_imp_per <- build_combined_data(c("perception", "explicit", "implicit"))
plot_combined(
  data_exp_imp_per,
  "Automaticity by pathway (explicit/implicit/perception)",
  out_exp_imp_per
)
