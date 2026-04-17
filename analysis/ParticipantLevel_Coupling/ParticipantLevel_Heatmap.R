library(dplyr)
library(tidyr)
library(ggplot2)

default_r <- "data/processed/participant_level_4x4_spearman_r.csv"
default_p <- "data/processed/participant_level_4x4_spearman_p.csv"
default_out <- "outputs/Trial_decoupling/participant_level_4x4_spearman_heatmap.png"

args <- commandArgs(trailingOnly = TRUE)
r_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_r
p_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_p
out_path <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_out

r_mat <- read.csv(r_path, row.names = 1, check.names = FALSE)
p_mat <- read.csv(p_path, row.names = 1, check.names = FALSE)

r_long <- r_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("row") %>%
  pivot_longer(-row, names_to = "col", values_to = "r")
p_long <- p_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("row") %>%
  pivot_longer(-row, names_to = "col", values_to = "p")

df <- r_long %>%
  left_join(p_long, by = c("row", "col"))

label_map <- c(
  dorsal_exp = "Trajectory–Instructed",
  ventral_exp = "Luminance–Instructed",
  dorsal_imp = "Trajectory–Passive",
  ventral_imp = "Luminance–Passive"
)
label_order <- unname(label_map[c("dorsal_exp", "ventral_exp", "dorsal_imp", "ventral_imp")])

df <- df %>%
  mutate(
    row_label = recode(row, !!!label_map),
    col_label = recode(col, !!!label_map),
    is_diag = row == col,
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = ifelse(is_diag, "", ifelse(sig == "", sprintf("%.2f", r), paste0(sprintf("%.2f", r), "\n", sig))),
    row_label = factor(row_label, levels = label_order),
    col_label = factor(col_label, levels = label_order)
  )

df <- df %>%
  mutate(r_plot = ifelse(is_diag, NA_real_, r))

p <- ggplot(df, aes(x = col_label, y = row_label, fill = r_plot)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3.5, color = "black") +
  scale_fill_gradient2(
    low = "#F6D6D6",
    mid = "white",
    high = "#2B6CB0",
    midpoint = 0,
    name = "r",
    na.value = "#F2F2F2"
  ) +
  scale_y_discrete(limits = rev(label_order)) +
  labs(
    title = "Participant-level coupling (Spearman r)",
    x = NULL,
    y = NULL
  ) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(family = "Times New Roman"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)
  )

ggsave(out_path, p, width = 6, height = 5.2, dpi = 300)
cat("Saved:", out_path, "\n")
