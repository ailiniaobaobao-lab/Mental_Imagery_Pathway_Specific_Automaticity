library(dplyr)
library(tidyr)

# Argument 1: coupling trial CSV. Argument 2: output directory.
default_input <- "outputs/Trial_decoupling/automaticity_trial_coupling_data.csv"
default_output_dir <- "outputs/Trial_decoupling"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_dir <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_csv <- file.path(output_dir, "mechanistic_link_coupling_automaticity.csv")
output_txt <- file.path(output_dir, "mechanistic_link_coupling_automaticity.txt")

# Load trial-level coupling data and keep the phases used in the mechanistic test.
df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = tolower(phase),
    stimulus = tolower(stimulus),
    participant = as.character(participant),
    trial_id = as.character(trial_id),
    trial_key = paste(stimulus, trial_id, sep = "_")
  ) %>%
  filter(phase %in% c("explicit", "implicit"))

calc_r <- function(x, y) {
  if (sum(is.finite(x) & is.finite(y)) < 2) return(NA_real_)
  suppressWarnings(stats::cor(x, y, use = "complete.obs", method = "spearman"))
}

# Residualize both pathways within each participant and stimulus before estimating coupling.
df_resid <- df %>%
  group_by(participant, stimulus) %>%
  mutate(
    corridor_c = corridor_z - mean(corridor_z, na.rm = TRUE),
    ventral_c = ventral_z - mean(ventral_z, na.rm = TRUE)
  ) %>%
  ungroup()

# Summarize coupling strength separately for explicit and implicit phases.
coupling_by_phase <- df_resid %>%
  group_by(participant, phase) %>%
  summarize(
    r_couple = calc_r(corridor_c, ventral_c),
    med_couple = median(corridor_c * ventral_c, na.rm = TRUE),
    n_trials = n_distinct(trial_key),
    .groups = "drop"
  )

coupling_wide <- coupling_by_phase %>%
  pivot_wider(
    names_from = phase,
    values_from = c(r_couple, med_couple, n_trials),
    names_sep = "_"
  ) %>%
  rename(
    r_couple_exp = r_couple_explicit,
    r_couple_imp = r_couple_implicit,
    med_couple_exp = med_couple_explicit,
    med_couple_imp = med_couple_implicit,
    n_trials_exp = n_trials_explicit,
    n_trials_imp = n_trials_implicit
  )

coupling_overall <- df_resid %>%
  group_by(participant) %>%
  summarize(
    r_couple_overall = calc_r(corridor_c, ventral_c),
    med_couple_overall = median(corridor_c * ventral_c, na.rm = TRUE),
    n_trials_overall = n_distinct(trial_key),
    .groups = "drop"
  )

# Convert both pathways to within-participant z scores, then derive phase-level automaticity indices.
df_z <- df %>%
  group_by(participant) %>%
  mutate(
    dorsal_z = (corridor_z - mean(corridor_z, na.rm = TRUE)) / sd(corridor_z, na.rm = TRUE),
    ventral_z = (ventral_z - mean(ventral_z, na.rm = TRUE)) / sd(ventral_z, na.rm = TRUE)
  ) %>%
  ungroup()

phase_means <- df_z %>%
  group_by(participant, phase) %>%
  summarize(
    dorsal = mean(dorsal_z, na.rm = TRUE),
    ventral = mean(ventral_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = phase, values_from = c(dorsal, ventral), names_sep = "_")

metrics <- phase_means %>%
  mutate(
    delta_dorsal = dorsal_explicit - dorsal_implicit,
    delta_ventral = ventral_explicit - ventral_implicit,
    delta_diff = delta_dorsal - delta_ventral
  ) %>%
  left_join(coupling_wide, by = "participant") %>%
  left_join(coupling_overall, by = "participant")

write.csv(metrics, output_csv, row.names = FALSE)

# Fit the main and sensitivity models linking coupling to automaticity.
metrics_use <- metrics %>%
  filter(is.finite(med_couple_imp), is.finite(delta_diff))

model_main <- lm(delta_diff ~ med_couple_imp, data = metrics_use)
model_ctrl <- lm(delta_diff ~ med_couple_imp + n_trials_imp, data = metrics_use)
model_dorsal <- lm(delta_dorsal ~ med_couple_imp, data = metrics_use)
model_ventral <- lm(delta_ventral ~ med_couple_imp, data = metrics_use)
model_robust <- lm(delta_diff ~ r_couple_imp, data = metrics_use)

summary_lines <- c(
  "=== Mechanistic link: coupling vs automaticity ===",
  "",
  "Main: delta_diff ~ med_couple_imp",
  capture.output(summary(model_main)),
  "",
  "Control: delta_diff ~ med_couple_imp + n_trials_imp",
  capture.output(summary(model_ctrl)),
  "",
  "Secondary: delta_dorsal ~ med_couple_imp",
  capture.output(summary(model_dorsal)),
  "",
  "Secondary: delta_ventral ~ med_couple_imp",
  capture.output(summary(model_ventral)),
  "",
  "Robustness: delta_diff ~ r_couple_imp",
  capture.output(summary(model_robust))
)

writeLines(summary_lines, con = output_txt)
cat("Saved:", output_csv, "\n")
cat("Saved:", output_txt, "\n")
