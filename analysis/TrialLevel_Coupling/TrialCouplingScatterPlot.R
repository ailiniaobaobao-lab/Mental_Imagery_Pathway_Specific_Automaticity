library(dplyr)
library(ggplot2)


default_input <- "data/processed/automaticity_trial_coupling_data.csv"
default_output <- "outputs/Trial_decoupling/automaticity_trial_coupling_scatter.png"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output

df <- read.csv(input_path, stringsAsFactors = FALSE) %>%
  mutate(
    phase = tolower(phase),
    corridor_z_w = as.numeric(corridor_z_w),
    ventral_z_w = as.numeric(ventral_z_w)
  ) %>%
  filter(phase %in% c("explicit", "implicit")) %>%
  filter(is.finite(corridor_z_w), is.finite(ventral_z_w))

df <- df %>%
  mutate(phase = factor(phase, levels = c("explicit", "implicit")))

limits <- c(-25, 25)

line_fit <- lm(ventral_z_w ~ corridor_z_w, data = df)
line_intercept <- unname(coef(line_fit)[1])
line_slope <- unname(coef(line_fit)[2])

rmcorr_label <- "rmcorr not available"
if (requireNamespace("rmcorr", quietly = TRUE)) {
  rmc <- try(rmcorr::rmcorr(participant, corridor_z, ventral_z, data = df), silent = TRUE)
  if (!inherits(rmc, "try-error")) {
    r_label <- format(rmc$r, digits = 3)
    p_label <- ifelse(rmc$p < 0.001, format(rmc$p, digits = 2, scientific = TRUE), format(rmc$p, digits = 3))
    rmcorr_label <- paste0("rmcorr: r = ", r_label, ", p = ", p_label)
  }
}

p <- ggplot(df, aes(x = corridor_z_w, y = ventral_z_w)) +
  geom_point(alpha = 0.6, size = 1.7, color = "grey30") +
  geom_abline(intercept = line_intercept, slope = line_slope, color = "#2E6F95", linewidth = 0.9) +
  coord_cartesian(ylim = limits) +
  labs(
    x = "Trajectory (centered hit rate)",
    y = "Luminance (centered ventral_z)",
    caption = "Trimmed axis for visualization (y limited to -25 to 25)"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.02,
    vjust = 1.2,
    size = 3,
    color = "grey20",
    label = rmcorr_label,
    family = "Times New Roman"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.caption = element_text(size = 9, hjust = 0),
    text = element_text(family = "Times New Roman")
  )

ggsave(output_path, p, width = 6, height = 4.6, dpi = 300)
cat("Saved:", output_path, "\n")
