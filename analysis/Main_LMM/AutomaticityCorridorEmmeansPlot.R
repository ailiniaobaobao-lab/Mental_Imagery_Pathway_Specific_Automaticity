library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(emmeans)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}

# Argument 1: ventral trial CSV. Argument 2: corridor trial CSV. Argument 3: optional output directory.
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_corridor <- "data/processed/gaze_corridor_hit_rate_trials.csv"
default_output_dir <- "outputs/Gaze_corrodor"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
corridor_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_corridor
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_png <- file.path(output_dir, "figure2_automaticity_corridor_emmeans.png")
output_csv <- file.path(output_dir, "figure2_automaticity_corridor_emmeans.csv")

exclude_participants <- c("11191401", "11191503", "11210905", "11241311", "11241210", "03041338")

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

exclude_participants <- pad_id(exclude_participants)

# ==== Load Data ====
ventral <- read.csv(ventral_path, stringsAsFactors = FALSE)
corridor <- read.csv(corridor_path, stringsAsFactors = FALSE)

# ==== Ventral (perception-z, signed) ====
ventral_df_all <- ventral %>%
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

ventral_signed_all <- ventral_df_all %>%
  transmute(
    participant,
    phase,
    stimulus,
    pathway = "ventral",
    dv_raw = ifelse(stimulus == "white", -pupil_response, pupil_response)
  )

ventral_perc_stats_z_signed <- ventral_signed_all %>%
  filter(phase == "perception") %>%
  group_by(participant, stimulus) %>%
  summarize(
    base_mean = mean(dv_raw, na.rm = TRUE),
    base_sd = sd(dv_raw, na.rm = TRUE),
    .groups = "drop"
  )

ventral_z_all_signed <- ventral_signed_all %>%
  left_join(ventral_perc_stats_z_signed, by = c("participant", "stimulus")) %>%
  mutate(dv_z = (dv_raw - base_mean) / base_sd) %>%
  filter(is.finite(dv_z))

# ==== Corridor (raw) ====
corridor_df_all <- corridor %>%
  mutate(
    phase = factor(phase, levels = c("explicit", "implicit", "perception")),
    stimulus = factor(stimulus, levels = c("white", "black")),
    participant = pad_id(participant)
  ) %>%
  filter(!is.na(hit_rate), !is.na(participant), !is.na(phase), !is.na(stimulus)) %>%
  filter(!participant %in% exclude_participants)

corridor_for_combined <- corridor_df_all %>%
  transmute(
    participant,
    phase,
    stimulus,
    pathway = "dorsal",
    dv_raw = hit_rate
  )

# Build the shared long-format dataset used for the emmeans comparison plot.
build_combined <- function(ventral_df, dorsal_df, phases) {
  ventral_use <- ventral_df %>%
    filter(phase %in% phases) %>%
    transmute(
      participant,
      phase,
      stimulus,
      pathway = "ventral",
      dv_z
    )
  dorsal_use <- dorsal_df %>%
    filter(phase %in% phases) %>%
    transmute(
      participant,
      phase,
      stimulus,
      pathway = "dorsal",
      dv_raw
    )
  combined <- bind_rows(ventral_use %>% rename(dv_raw = dv_z), dorsal_use) %>%
    group_by(participant, pathway) %>%
    mutate(
      dv_mean = mean(dv_raw, na.rm = TRUE),
      dv_sd = sd(dv_raw, na.rm = TRUE),
      dv_z = (dv_raw - dv_mean) / dv_sd
    ) %>%
    ungroup() %>%
    filter(is.finite(dv_z)) %>%
    mutate(
      participant = factor(participant),
      pathway = factor(pathway, levels = c("ventral", "dorsal")),
      phase = factor(phase, levels = phases),
      stimulus = factor(stimulus, levels = c("white", "black"))
    )
  combined
}

combined_signed_corridor <- build_combined(ventral_z_all_signed, corridor_for_combined, c("explicit", "implicit"))

model <- lmer(dv_z ~ phase * pathway + stimulus + (1 | participant), data = combined_signed_corridor, REML = FALSE)

emm <- emmeans(model, ~ phase * pathway)
emm_df <- as.data.frame(emm)
write.csv(emm_df, output_csv, row.names = FALSE)

emm_df <- emm_df %>%
  mutate(
    phase = factor(phase, levels = c("explicit", "implicit")),
    pathway = factor(pathway, levels = c("ventral", "dorsal"))
  )

coef_df <- summary(model)$coefficients
beta <- NA_real_
p_val <- NA_real_
term_name <- "phaseimplicit:pathwaydorsal"
if (term_name %in% rownames(coef_df)) {
  beta <- coef_df[term_name, "Estimate"]
  p_val <- coef_df[term_name, "Pr(>|t|)"]
}
annot_text <- if (is.finite(beta) && is.finite(p_val)) {
  sprintf("Phase × Signature: \u03b2 = %.3f, p = %.3f", beta, p_val)
} else {
  "Phase \u00d7 Signature: \u03b2 = NA, p = NA"
}

phase_labels <- c(
  explicit = "Instructed Imagery",
  implicit = "Passive Viewing"
)

pathway_labels <- c(
  ventral = "Luminance signature",
  dorsal = "Trajectory signature"
)


p <- ggplot(emm_df, aes(x = phase, y = emmean, color = pathway, group = pathway)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.08, linewidth = 0.6) +
  labs(
    x = NULL,
    y = "Standardized signature strength (within-subject z)",
    color = "Signature"
  ) +
  scale_x_discrete(labels = phase_labels) +
  scale_color_manual(values = c(ventral = "#1F77B4", dorsal = "#D62728"), labels = pathway_labels) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = annot_text,
    hjust = 1.02,
    vjust = 1.1,
    size = 3,
    color = "grey20",
    family = "Times New Roman"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    text = element_text(family = "Times New Roman")
  )

ggsave(output_png, p, width = 7.5, height = 5.2, dpi = 300)
cat("Saved:", output_png, "\n")
cat("Saved:", output_csv, "\n")
