library(dplyr)
library(ggplot2)
library(readxl)

# input
default_ventral <- "data/processed/ventral_trial_responses.csv"
default_info <- "data/metadata/Participants_info.xlsx"
default_output <- "outputs/PupilSize_Raw/ventral_phase_group_raincloud.png"

args <- commandArgs(trailingOnly = TRUE)
ventral_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_ventral
info_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_info
output_path <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output

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

add_perception_z <- function(df, dv_col) {
  df <- df %>%
    filter(phase %in% c("perception", "explicit", "implicit")) %>%
    filter(!is.na(.data[[dv_col]]), !is.na(participant), !is.na(stimulus))
  base <- df %>%
    filter(phase == "perception") %>%
    group_by(participant, stimulus) %>%
    summarize(
      base_mean = mean(.data[[dv_col]], na.rm = TRUE),
      base_sd = sd(.data[[dv_col]], na.rm = TRUE),
      .groups = "drop"
    )
  df %>%
    left_join(base, by = c("participant", "stimulus")) %>%
    mutate(dv_z = (.data[[dv_col]] - base_mean) / base_sd) %>%
    filter(is.finite(dv_z))
}

info <- readxl::read_excel(info_path)
col_h <- names(info)[8]
group_map <- info %>%
  transmute(
    participant = pad_id(.data[["subject#"]]),
    group = as.character(.data[[col_h]])
  ) %>%
  filter(group %in% c("B", "E", "N")) %>%
  mutate(group = ifelse(group == "E", "I", group))

ventral <- read.csv(ventral_path, stringsAsFactors = FALSE) %>%
  mutate(
    participant = pad_id(participant),
    phase = normalize_phase(phase),
    stimulus = tolower(stimulus),
    dv_raw = ifelse(stimulus == "white", -pupil_response, pupil_response)
  )

ventral_z <- add_perception_z(ventral, "dv_raw")

ventral_sum <- ventral_z %>%
  group_by(participant, phase) %>%
  summarize(ventral_value = mean(dv_z, na.rm = TRUE), .groups = "drop") %>%
  inner_join(group_map, by = "participant") %>%
  filter(participant != "02131225")

phase_labels <- c(
  perception = "Perception",
  explicit = "Instructed Imagery",
  implicit = "Passive Viewing"
)

ventral_sum <- ventral_sum %>%
  mutate(
    phase = factor(phase, levels = c("perception", "explicit", "implicit")),
    group = factor(group, levels = c("B", "I", "N"))
  )

pd <- position_dodge(width = 0.75)
pj <- position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75)

p <- ggplot(ventral_sum, aes(x = phase, y = ventral_value, color = group, fill = group)) +
  {if (requireNamespace("gghalves", quietly = TRUE)) {
    gghalves::geom_half_violin(
      side = "l",
      position = pd,
      alpha = 0.25,
      trim = FALSE
    )
  } else {
    geom_violin(
      position = pd,
      alpha = 0.25,
      trim = FALSE
    )
  }} +
  geom_boxplot(width = 0.12, position = pd, outlier.shape = NA, alpha = 0.5) +
  geom_point(position = pj, size = 1.4, alpha = 0.7) +
  scale_x_discrete(labels = phase_labels) +
  labs(
    x = "Phase",
    y = "Luminance signature (perception-z)",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(text = element_text(family = "Times New Roman"))

ggsave(output_path, p, width = 7.2, height = 4.8, dpi = 300)
cat("Saved:", output_path, "\n")
