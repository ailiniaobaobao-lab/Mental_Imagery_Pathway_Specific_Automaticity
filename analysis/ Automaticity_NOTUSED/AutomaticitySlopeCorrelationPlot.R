library(dplyr)
library(ggplot2)

# ==== Input/Output ====
# Argument 1: slopes CSV. Argument 2: optional output directory. Argument 3: optional bootstrap count.
# Argument 4: rope_r (default 0.3). Argument 5: alpha (default 0.05).
default_slopes <- "data/processed/automaticity_onoff_slopes.csv"
default_output_dir <- "outputs/PupilSize_Raw"

args <- commandArgs(trailingOnly = TRUE)
slopes_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_slopes
output_dir <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output_dir
boot_n <- if (length(args) >= 3) as.integer(args[3]) else 2000L
rope_r <- if (length(args) >= 4) as.numeric(args[4]) else 0.30
alpha <- if (length(args) >= 5) as.numeric(args[5]) else 0.05

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_plot <- file.path(output_dir, "automaticity_slope_correlation.png")
output_stats <- file.path(output_dir, "automaticity_slope_correlation.csv")
plot_x_max <- 0.35

slopes <- read.csv(slopes_path, stringsAsFactors = FALSE) %>%
  filter(!is.na(slope_ventral), !is.na(slope_dorsal))

if (!nrow(slopes)) {
  stop("No valid slope rows in: ", slopes_path)
}

set.seed(202401)
r_obs <- cor(slopes$slope_ventral, slopes$slope_dorsal, use = "complete.obs")

boot_r <- replicate(boot_n, {
  idx <- sample.int(nrow(slopes), replace = TRUE)
  cor(slopes$slope_ventral[idx], slopes$slope_dorsal[idx], use = "complete.obs")
})

ci <- stats::quantile(boot_r, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 7)

fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
z_obs <- fisher_z(r_obs)
z_low <- fisher_z(-abs(rope_r))
z_high <- fisher_z(abs(rope_r))
se <- 1 / sqrt(nrow(slopes) - 3)
z1 <- (z_obs - z_low) / se
z2 <- (z_obs - z_high) / se
ptost_lower <- 1 - pnorm(z1)
ptost_upper <- pnorm(z2)
equivalent <- (ptost_lower < alpha) && (ptost_upper < alpha)

stats_df <- data.frame(
  n = nrow(slopes),
  r = r_obs,
  ci_low = ci[1],
  ci_high = ci[2],
  boot_n = boot_n,
  rope_r = rope_r,
  tost_p_lower = ptost_lower,
  tost_p_upper = ptost_upper,
  equivalent = equivalent,
  alpha = alpha
)

subtitle_text <- sprintf("r = %.3f  (95%% CI %.3f to %.3f), n = %d", r_obs, ci[1], ci[2], nrow(slopes))
out_of_range <- sum(slopes$slope_ventral > plot_x_max)
out_label <- if (out_of_range > 0) "One point falls outside the x-axis range" else NULL

p <- ggplot(slopes, aes(x = slope_ventral, y = slope_dorsal)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#2C7FB8") +
  coord_cartesian(xlim = c(NA, plot_x_max)) +
  labs(
    title = "Ventral vs Dorsal Phase Slopes",
    subtitle = subtitle_text,
    x = "Slope (ventral: implicit vs explicit)",
    y = "Slope (dorsal: implicit vs explicit)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

if (!is.null(out_label)) {
  p <- p + annotate("text",
    x = plot_x_max,
    y = Inf,
    label = out_label,
    hjust = 1.02,
    vjust = 1.2,
    size = 3.5,
    color = "grey30"
  )
}

ggsave(output_plot, p, width = 6.5, height = 5, dpi = 300)
write.csv(stats_df, output_stats, row.names = FALSE)

cat("Saved:", output_plot, "\n")
cat("Saved:", output_stats, "\n")
