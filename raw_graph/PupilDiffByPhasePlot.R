library(readxl)
library(dplyr)
library(ggplot2)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}

# Argument 1: ventral CSV. Argument 2: participants xlsx. Argument 3: output path.
default_input <- "data/processed/ventral_trial_responses.csv"
default_info <- "data/metadata/Participants_info.xlsx"
default_output <- "outputs/PupilSize_Raw/pupil_response_diff_bw_by_phase.png"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_path <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output

normalize_phase <- function(x) {
  dplyr::case_when(
    x == "phase_p" ~ "perception",
    x == "phase_e" ~ "explicit",
    x == "phase_i" ~ "implicit",
    TRUE ~ as.character(x)
  )
}

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

# Load the checked participant list so the phase trend only uses included cases.
info <- readxl::read_excel(info_path)
checked_ids <- info %>%
  filter(.data[["...9"]] == TRUE) %>%
  transmute(participant = pad_id(.data[["subject#"]])) %>%
  pull(participant)

# Load ventral trial summaries and compute participant-level black-minus-white differences by phase.
df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = normalize_phase(phase)
  )

df_sum <- df %>%
  filter(participant %in% checked_ids) %>%
  group_by(participant, phase) %>%
  summarize(pupil_response_diff_bw = mean(pupil_response_diff_bw, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(pupil_response_diff_bw))

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

# Test the linear phase trend with a mixed model.
model <- lmer(pupil_response_diff_bw ~ phase_num + (1 | participant), data = df_sum, REML = FALSE)
coef_df <- summary(model)$coefficients
if ("phase_num" %in% rownames(coef_df) && "Pr(>|t|)" %in% colnames(coef_df)) {
  p_val <- coef_df["phase_num", "Pr(>|t|)"]
} else {
  p_val <- NA_real_
}

# Plot the participant trajectories and annotate the phase-trend result.
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

p <- ggplot(df_sum, aes(x = phase, y = pupil_response_diff_bw, group = participant)) +
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
    title = "Pupil response diff (black - white) by phase",
    x = "Phase",
    y = "pupil_response_diff_bw"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(output_path, p, width = 6.5, height = 4.8, dpi = 300)
cat("Saved:", output_path, "\n")
