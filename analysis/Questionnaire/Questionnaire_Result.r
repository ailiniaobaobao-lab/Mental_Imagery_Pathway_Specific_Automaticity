library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# input path
default_input <- "data/metadata/Participants_info.xlsx"
default_output_dir <- "outputs"

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_dir <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_csv <- file.path(output_dir, "Participants_info_correlations.csv")
output_png <- file.path(output_dir, "Participants_info_corr_heatmap.png")

cols <- c("VVIQ", "VOSI_O", "VOSI_S", "VOSI_total")
df <- readxl::read_excel(input_path)
sub <- df %>% select(all_of(cols))

# correlation chart
pairs <- combn(cols, 2, simplify = FALSE)
out <- lapply(pairs, function(p) {
  x <- sub[[p[1]]]
  y <- sub[[p[2]]]
  keep <- complete.cases(x, y)
  if (sum(keep) < 3) {
    return(data.frame(var1 = p[1], var2 = p[2], n = sum(keep), r = NA, p = NA))
  }
  ct <- cor.test(x[keep], y[keep], method = "pearson")
  data.frame(
    var1 = p[1],
    var2 = p[2],
    n = sum(keep),
    r = unname(ct[["estimate"]]),
    p = ct[["p.value"]]
  )
})
res <- bind_rows(out)
write.csv(res, output_csv, row.names = FALSE)

# heatmap
cor_mat <- cor(sub, use = "pairwise.complete.obs", method = "pearson")
cor_df <- as.data.frame(cor_mat)
cor_df$var1 <- rownames(cor_df)
long <- pivot_longer(cor_df, -var1, names_to = "var2", values_to = "r")
long$var1 <- factor(long$var1, levels = cols)
long$var2 <- factor(long$var2, levels = cols)

p <- ggplot(long, aes(x = var1, y = var2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#2C7FB8",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_equal() +
  labs(
    title = "Correlation Heatmap (Pearson r)",
    x = NULL,
    y = NULL,
    fill = "r"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))

ggsave(output_png, p, width = 5.5, height = 5, dpi = 300)

cat("saved:", output_csv, "\n")
cat("saved:", output_png, "\n")
