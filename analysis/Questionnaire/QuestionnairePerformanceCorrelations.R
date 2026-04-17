library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)

default_info <- "data/metadata/Participants_info.xlsx"
default_perf <- "outputs/Trial_decoupling/within_phase_participant_coupling_summary.csv"
default_output_dir <- "outputs/Questionnaire_analysis"

args <- commandArgs(trailingOnly = TRUE)
info_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_info
perf_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = TRUE) else default_perf
output_dir <- if (length(args) >= 3) normalizePath(args[3], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_csv <- file.path(output_dir, "questionnaire_performance_correlations.csv")
output_txt <- file.path(output_dir, "questionnaire_performance_correlations.txt")
output_plot <- file.path(output_dir, "questionnaire_performance_forest.png")
output_plot_vosi_di <- file.path(output_dir, "questionnaire_vosi_spatial_vs_di.png")
output_plot_vviq_vi <- file.path(output_dir, "questionnaire_vviq_vs_vi.png")
output_plot_vosi_o_vi <- file.path(output_dir, "questionnaire_vosi_o_vs_vi.png")
output_plot_vviq_ve <- file.path(output_dir, "questionnaire_vviq_vs_ve.png")
output_plot_vosi_o_ve <- file.path(output_dir, "questionnaire_vosi_o_vs_ve.png")
output_plot_vosi_s_ve <- file.path(output_dir, "questionnaire_vosi_s_vs_ve.png")
output_plot_vviq_de <- file.path(output_dir, "questionnaire_vviq_vs_de.png")
output_plot_vosi_o_de <- file.path(output_dir, "questionnaire_vosi_o_vs_de.png")
output_plot_vosi_s_de <- file.path(output_dir, "questionnaire_vosi_s_vs_de.png")

# Optional per-plot participant exclusions; keep empty unless a specific rerun requires them.
exclude_vi_participants <- character(0)
exclude_ve_participants <- character(0)

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

info <- readxl::read_excel(info_path)
subj_col <- names(info)[tolower(names(info)) %in% c("subject#", "subject", "participant")]
if (!length(subj_col)) stop("No subject column found in Participants_info.xlsx")
subj_col <- subj_col[1]

if (ncol(info) < 6) stop("Participants_info.xlsx has fewer than 6 columns; cannot read D/E/F.")
q_cols <- names(info)[4:6]

q_df <- info %>%
  transmute(
    participant = pad_id(.data[[subj_col]]),
    VVIQ_total = suppressWarnings(as.numeric(.data[[q_cols[1]]])),
    VOSI_O = suppressWarnings(as.numeric(.data[[q_cols[2]]])),
    VOSI_S = suppressWarnings(as.numeric(.data[[q_cols[3]]]))
  )

perf_long <- read.csv(perf_path, stringsAsFactors = FALSE)
perf_wide <- perf_long %>%
  mutate(
    participant = pad_id(participant),
    phase = tolower(phase)
  ) %>%
  filter(phase %in% c("explicit", "implicit")) %>%
  pivot_wider(
    names_from = phase,
    values_from = c(dorsal_mean, ventral_mean),
    names_sep = "_"
  ) %>%
  transmute(
    participant,
    DE = dorsal_mean_explicit,
    DI = dorsal_mean_implicit,
    VE = ventral_mean_explicit,
    VI = ventral_mean_implicit
  )

df <- perf_wide %>%
  inner_join(q_df, by = "participant")

cor_row <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < 3) return(list(r = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(x, y, method = method))
  list(r = unname(ct$estimate), p = ct$p.value, n = n)
}

questionnaires <- c("VVIQ_total", "VOSI_O", "VOSI_S")
performances <- c("DE", "DI", "VE", "VI")

results <- bind_rows(lapply(questionnaires, function(q) {
  bind_rows(lapply(performances, function(p) {
    pear <- cor_row(df[[q]], df[[p]], "pearson")
    spear <- cor_row(df[[q]], df[[p]], "spearman")
    tibble(
      questionnaire = q,
      performance = p,
      pearson_r = pear$r,
      pearson_p = pear$p,
      pearson_n = pear$n,
      spearman_r = spear$r,
      spearman_p = spear$p,
      spearman_n = spear$n
    )
  }))
}))

write.csv(results, output_csv, row.names = FALSE)

fmt <- function(x, digits = 3) {
  if (!is.finite(x)) return("NA")
  formatC(x, digits = digits, format = "f")
}

fisher_ci <- function(r, n, level = 0.95) {
  if (!is.finite(r) || n < 4 || abs(r) >= 1) {
    return(c(NA_real_, NA_real_))
  }
  z <- atanh(r)
  se <- 1 / sqrt(n - 3)
  z_crit <- qnorm(1 - (1 - level) / 2)
  ci_low <- tanh(z - z_crit * se)
  ci_high <- tanh(z + z_crit * se)
  c(ci_low, ci_high)
}

plot_df <- results %>%
  mutate(
    questionnaire_label = recode(
      questionnaire,
      VVIQ_total = "VVIQ",
      VOSI_O = "VOSI-O",
      VOSI_S = "VOSI-S"
    ),
    label = paste0(questionnaire_label, "-", performance),
    r = pearson_r,
    n = pearson_n,
    ci = purrr::map2(r, n, fisher_ci),
    ci_low = purrr::map_dbl(ci, 1),
    ci_high = purrr::map_dbl(ci, 2)
  ) %>%
  arrange(questionnaire_label, performance) %>%
  mutate(label = factor(label, levels = rev(unique(label))))

p <- ggplot(plot_df, aes(x = r, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, color = "gray40") +
  geom_point(size = 2, color = "black") +
  labs(x = "Pearson r", y = "Questionnaire × Index") +
  theme_classic(base_size = 12)

ggsave(output_plot, p, width = 8, height = 6, dpi = 300)

make_scatter <- function(df_in, x_col, y_col, x_label, y_label, out_path, xlim = NULL, ylim = NULL) {
  keep <- is.finite(df_in[[x_col]]) & is.finite(df_in[[y_col]])
  df_plot <- df_in[keep, , drop = FALSE]
  ct <- suppressWarnings(cor.test(df_plot[[x_col]], df_plot[[y_col]], method = "pearson"))
  ann <- paste0("r = ", fmt(unname(ct$estimate)), ", p = ", fmt(ct$p.value))
  g <- ggplot(df_plot, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(color = "black", size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "gray40", fill = "gray80") +
    annotate("text", x = Inf, y = Inf, label = ann, hjust = 1.05, vjust = 1.2, size = 3.5) +
    labs(x = x_label, y = y_label) +
    theme_classic(base_size = 12) +
    theme(text = element_text(family = "Times New Roman"))
  if (!is.null(xlim) || !is.null(ylim)) {
    g <- g + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  ggsave(out_path, g, width = 6, height = 4.5, dpi = 300)
}

df_vi <- df %>% filter(!participant %in% exclude_vi_participants)
df_ve <- df %>% filter(!participant %in% exclude_ve_participants)
ventral_ylim <- range(c(df_vi$VI, df_ve$VE), na.rm = TRUE)
dorsal_ylim <- range(c(df$DE, df$DI), na.rm = TRUE)
questionnaire_xlim <- c(1, 5)
make_scatter(
  df,
  "VOSI_S",
  "DI",
  "VOSI-Spatial",
  "Trajectory–Passive",
  output_plot_vosi_di,
  xlim = questionnaire_xlim,
  ylim = dorsal_ylim
)
make_scatter(
  df_vi,
  "VVIQ_total",
  "VI",
  "VVIQ",
  "Luminance–Passive",
  output_plot_vviq_vi,
  xlim = questionnaire_xlim,
  ylim = ventral_ylim
)
make_scatter(
  df_vi,
  "VOSI_O",
  "VI",
  "VOSI-Object",
  "Luminance–Passive",
  output_plot_vosi_o_vi,
  xlim = questionnaire_xlim,
  ylim = ventral_ylim
)
make_scatter(
  df_ve,
  "VVIQ_total",
  "VE",
  "VVIQ",
  "Luminance–Instructed",
  output_plot_vviq_ve,
  xlim = questionnaire_xlim,
  ylim = ventral_ylim
)
make_scatter(
  df_ve,
  "VOSI_O",
  "VE",
  "VOSI-Object",
  "Luminance–Instructed",
  output_plot_vosi_o_ve,
  xlim = questionnaire_xlim,
  ylim = ventral_ylim
)
make_scatter(
  df_ve,
  "VOSI_S",
  "VE",
  "VOSI-Spatial",
  "Luminance–Instructed",
  output_plot_vosi_s_ve,
  xlim = questionnaire_xlim,
  ylim = ventral_ylim
)
make_scatter(
  df,
  "VVIQ_total",
  "DE",
  "VVIQ",
  "Trajectory–Instructed",
  output_plot_vviq_de,
  xlim = questionnaire_xlim,
  ylim = dorsal_ylim
)
make_scatter(
  df,
  "VOSI_O",
  "DE",
  "VOSI-Object",
  "Trajectory–Instructed",
  output_plot_vosi_o_de,
  xlim = questionnaire_xlim,
  ylim = dorsal_ylim
)
make_scatter(
  df,
  "VOSI_S",
  "DE",
  "VOSI-Spatial",
  "Trajectory–Instructed",
  output_plot_vosi_s_de,
  xlim = questionnaire_xlim,
  ylim = dorsal_ylim
)

summary_lines <- c(
  "=== Questionnaire vs performance correlations ===",
  paste0("Questionnaire columns (D/E/F): ", paste(q_cols, collapse = ", ")),
  ""
)

for (q in questionnaires) {
  summary_lines <- c(summary_lines, paste0("[", q, "]"))
  for (p in performances) {
    row <- results %>% filter(questionnaire == q, performance == p)
    summary_lines <- c(
      summary_lines,
      paste0(
        "  ", p,
        " | Pearson r=", fmt(row$pearson_r), ", p=", fmt(row$pearson_p), ", n=", row$pearson_n,
        " | Spearman r=", fmt(row$spearman_r), ", p=", fmt(row$spearman_p), ", n=", row$spearman_n
      )
    )
  }
  summary_lines <- c(summary_lines, "")
}

writeLines(summary_lines, con = output_txt)
cat("Saved:", output_csv, "\n")
cat("Saved:", output_txt, "\n")
cat("Saved:", output_plot, "\n")
cat("Saved:", output_plot_vosi_di, "\n")
cat("Saved:", output_plot_vviq_vi, "\n")
cat("Saved:", output_plot_vosi_o_vi, "\n")
cat("Saved:", output_plot_vviq_ve, "\n")
cat("Saved:", output_plot_vosi_o_ve, "\n")
cat("Saved:", output_plot_vosi_s_ve, "\n")
cat("Saved:", output_plot_vviq_de, "\n")
cat("Saved:", output_plot_vosi_o_de, "\n")
cat("Saved:", output_plot_vosi_s_de, "\n")
