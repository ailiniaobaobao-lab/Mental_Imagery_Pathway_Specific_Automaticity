library(readxl)
library(dplyr)
library(ggplot2)
library(lme4)
if (requireNamespace("lmerTest", quietly = TRUE)) {
  library(lmerTest)
}

# Argument 1: dorsal trial CSV. Argument 2: participants xlsx. Argument 3: output directory.
default_input <- "data/processed/dorsal_trial_errors.csv"
default_info <- "data/metadata/Participants_info.xlsx"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

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

# Load participant group labels used to stratify the dorsal phase trend plots.
info <- readxl::read_excel(info_path)
col_h <- names(info)[8]
group_map <- info %>%
  transmute(
    participant = pad_id(.data[["subject#"]]),
    group = as.character(.data[[col_h]])
  ) %>%
  filter(group %in% c("B", "N", "E"))

# Load dorsal trial data and average RMSE within participant and phase.
df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = normalize_phase(phase)
  )

df_sum <- df %>%
  group_by(participant, phase) %>%
  summarize(traj_rmse = mean(traj_rmse, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(traj_rmse)) %>%
  inner_join(group_map, by = "participant")

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

# Plot one phase trend figure for a single participant group.
plot_one <- function(group_label) {
  sub <- df_sum %>% filter(group == group_label)
  if (!nrow(sub)) return(NULL)
  p_val <- NA_real_
  model <- try(lmer(traj_rmse ~ phase_num + (1 | participant), data = sub, REML = FALSE), silent = TRUE)
  if (!inherits(model, "try-error")) {
    coef_df <- summary(model)$coefficients
    if ("phase_num" %in% rownames(coef_df) && "Pr(>|t|)" %in% colnames(coef_df)) {
      p_val <- coef_df["phase_num", "Pr(>|t|)"]
    }
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
  p <- ggplot(sub, aes(x = phase, y = traj_rmse, group = participant)) +
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
      title = paste0("Dorsal RMSE by phase (Group ", group_label, ")"),
      x = "Phase",
      y = "traj_rmse"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  out_path <- file.path(output_dir, paste0("dorsal_rmse_by_phase_group_", group_label, ".png"))
  ggsave(out_path, p, width = 6.5, height = 4.8, dpi = 300)
  cat("Saved:", out_path, "\n")
}

plot_one("B")
plot_one("N")
plot_one("E")
