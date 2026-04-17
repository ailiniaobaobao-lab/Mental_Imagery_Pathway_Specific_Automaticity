library(dplyr)
library(lme4)
library(emmeans)
library(readxl)
library(ggplot2)

# Argument 1: dorsal trial CSV. Argument 2: participants xlsx. Argument 3: output directory.
default_dorsal <- "data/processed/dorsal_trial_errors.csv"
default_info <- "data/metadata/Participants_info.xlsx"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
dorsal_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_dorsal
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_csv <- file.path(output_dir, "dorsal_group_emmeans.csv")
output_png <- file.path(output_dir, "dorsal_group_emmeans.png")

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

# Load participant group labels so dorsal emmeans can be estimated by subgroup.
info <- readxl::read_excel(info_path)
col_h <- names(info)[8]
group_map <- info %>%
  transmute(
    participant = pad_id(.data[["subject#"]]),
    group = as.character(.data[[col_h]])
  ) %>%
  filter(group %in% c("B", "N", "E"))

# Load dorsal trial data and convert RMSE into a strength-style metric.
strength_eps <- 1e-6
dorsal_raw <- read.csv(dorsal_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = normalize_phase(phase),
    stimulus = factor(stimulus, levels = c("white", "black")),
    strength_rmse = 1 / (pmax(traj_rmse, 0) + strength_eps)
  ) %>%
  filter(!is.na(strength_rmse), !is.na(participant), !is.na(phase), !is.na(stimulus)) %>%
  inner_join(group_map, by = "participant")

# Use perception trials to standardize each participant and stimulus combination.
perc_stats <- dorsal_raw %>%
  filter(phase == "perception") %>%
  group_by(participant, stimulus) %>%
  summarize(
    base_mean = mean(strength_rmse, na.rm = TRUE),
    base_sd = sd(strength_rmse, na.rm = TRUE),
    .groups = "drop"
  )

# Build the final analysis table and fit the mixed model used for marginal means.
dorsal <- dorsal_raw %>%
  left_join(perc_stats, by = c("participant", "stimulus")) %>%
  mutate(dv_z = (strength_rmse - base_mean) / base_sd) %>%
  filter(is.finite(dv_z)) %>%
  mutate(
    participant = factor(participant),
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    group = factor(group, levels = c("B", "N", "E"))
  )

model <- lmer(dv_z ~ phase * group + stimulus + (1 | participant), data = dorsal, REML = FALSE)

# Estimate marginal means and planned contrasts for the phase-by-group pattern.
emm_phase_group <- emmeans(model, ~ phase * group)
emm_df <- as.data.frame(emm_phase_group) %>%
  mutate(section = "emmeans")

emm_phase_by_group <- emmeans(model, ~ phase | group)
contrast_within <- contrast(
  emm_phase_by_group,
  method = list("implicit_vs_explicit" = c(-1, 1, 0)),
  by = "group",
  adjust = "none"
)
contrast_within_df <- as.data.frame(contrast_within) %>%
  mutate(section = "within_group_implicit_vs_explicit")

emm_group_implicit <- emmeans(model, ~ group, at = list(phase = "implicit"))
contrast_between <- contrast(emm_group_implicit, method = "pairwise", adjust = "none")
contrast_between_df <- as.data.frame(contrast_between) %>%
  mutate(section = "between_groups_at_implicit")

out_table <- bind_rows(
  emm_df,
  contrast_within_df,
  contrast_between_df
)

# Save the emmeans table and render the summary figure.
write.csv(out_table, output_csv, row.names = FALSE)

plot_df <- emm_df %>%
  mutate(
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    group = factor(group, levels = c("B", "N", "E"))
  )

p <- ggplot(plot_df, aes(x = phase, y = emmean, color = group, group = group)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, alpha = 0.7) +
  labs(
    title = "Dorsal group × phase (perception-z) estimated means",
    x = "Phase",
    y = "Estimated mean (dv_z)",
    color = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(output_png, p, width = 6.5, height = 4.8, dpi = 300)

cat("Saved:", output_csv, "\n")
cat("Saved:", output_png, "\n")
