library(dplyr)
library(ggplot2)
library(readxl)
library(scales)
library(tidyr)
library(grid)

# Argument 1: RDS directory. Argument 2: optional output file path.
default_rds_dir <- "data/processed/RDS_Data"
default_output <- "outputs/PupilSize_Raw/GazePlots_phase1_no_swap/gaze_phase1_white_combined.png"

participants <- c("02131225", "02021516", "02251433")
smooth_gaze <- TRUE
smooth_window <- 5
min_time_ms <- 1000
caption_text <- "Perception / Instructed Imagery / Passive Viewing"

args <- commandArgs(trailingOnly = TRUE)
rds_dir <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_rds_dir
output_path <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output
out_dir <- dirname(output_path)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==== Screen size ====
screen_width <- 4480
screen_height <- 2520
flip_coin_y <- TRUE
# ==== Content size ====
content_width <- 3840
content_height <- 2160
content_offset_x <- 0
content_offset_y <- 0

# ==== Coin paths (no swap) ====
coin_path_white_file <- "data/reference/White_Trajectory.xlsx"
coin_path_screen_w <- screen_width
coin_path_screen_h <- screen_height

occluder_white_file <- "data/reference/white_1_occluder.csv"
white_occluder_dx <- 0
white_occluder_dy <- -200

load_coin_path <- function(path, screen_width, screen_height, base_w, base_h, flip_y,
                           content_width, content_height, offset_x, offset_y) {
  if (!file.exists(path)) return(NULL)
  df <- readxl::read_excel(path)
  has_screen <- all(c("screen_x", "screen_y") %in% names(df))
  has_viewport <- all(c("viewport_x", "viewport_y") %in% names(df))
  if (!has_screen && !has_viewport) {
    warning("Coin path missing screen_x/screen_y or viewport_x/viewport_y: ", path)
    return(NULL)
  }
  if ("t" %in% names(df)) df <- df %>% arrange(t)
  if (has_screen) {
    base_w <- as.numeric(base_w)
    base_h <- as.numeric(base_h)
    if (is.na(base_w) || base_w <= 0) base_w <- screen_width
    if (is.na(base_h) || base_h <= 0) base_h <- screen_height
    sx <- screen_width / base_w
    sy <- screen_height / base_h
    out <- df %>%
      transmute(
        screen_x = as.numeric(screen_x) * sx,
        screen_y = as.numeric(screen_y) * sy
      ) %>%
      drop_na(screen_x, screen_y)
  } else {
    vw <- ifelse(is.finite(content_width) && content_width > 0, content_width, screen_width)
    vh <- ifelse(is.finite(content_height) && content_height > 0, content_height, screen_height)
    ox <- ifelse(is.finite(offset_x), offset_x, 0)
    oy <- ifelse(is.finite(offset_y), offset_y, 0)
    out <- df %>%
      transmute(
        screen_x = as.numeric(viewport_x) * vw + ox,
        screen_y = as.numeric(viewport_y) * vh + oy
      ) %>%
      drop_na(screen_x, screen_y)
  }
  if (isTRUE(flip_y)) {
    if (has_screen) {
      out <- out %>% mutate(screen_y = screen_height - screen_y)
    } else {
      out <- out %>% mutate(screen_y = (content_height - (screen_y - offset_y)) + offset_y)
    }
  }
  out
}

load_occluder <- function(path, screen_width, screen_height) {
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("left", "right", "bottom", "top") %in% names(df))) {
    warning("Occluder missing left/right/bottom/top: ", path)
    return(NULL)
  }
  row <- df[1, , drop = FALSE]
  screen_w <- if ("screen_w" %in% names(row)) as.numeric(row$screen_w) else screen_width
  screen_h <- if ("screen_h" %in% names(row)) as.numeric(row$screen_h) else screen_height
  if (is.na(screen_w) || screen_w == 0) screen_w <- screen_width
  if (is.na(screen_h) || screen_h == 0) screen_h <- screen_height
  sx <- screen_width / screen_w
  sy <- screen_height / screen_h
  left <- as.numeric(row$left) * sx
  right <- as.numeric(row$right) * sx
  bottom <- as.numeric(row$bottom) * sy
  top <- as.numeric(row$top) * sy
  bottom_t <- screen_height - bottom
  top_t <- screen_height - top
  xmin <- min(left, right, na.rm = TRUE)
  xmax <- max(left, right, na.rm = TRUE)
  ymin <- min(bottom_t, top_t, na.rm = TRUE)
  ymax <- max(bottom_t, top_t, na.rm = TRUE)
  xmin <- max(0, xmin)
  xmax <- min(screen_width, xmax)
  ymin <- max(0, ymin)
  ymax <- min(screen_height, ymax)
  if (is.na(xmin) || is.na(xmax) || is.na(ymin) || is.na(ymax) || xmin >= xmax || ymin >= ymax) {
    return(NULL)
  }
  data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

shift_occluder <- function(rect, dx, dy, screen_width, screen_height) {
  if (is.null(rect) || !nrow(rect)) return(NULL)
  shifted <- data.frame(
    xmin = rect$xmin + dx,
    xmax = rect$xmax + dx,
    ymin = rect$ymin + dy,
    ymax = rect$ymax + dy
  )
  shifted <- data.frame(
    xmin = max(0, min(shifted$xmin, shifted$xmax, na.rm = TRUE)),
    xmax = min(screen_width, max(shifted$xmin, shifted$xmax, na.rm = TRUE)),
    ymin = max(0, min(shifted$ymin, shifted$ymax, na.rm = TRUE)),
    ymax = min(screen_height, max(shifted$ymin, shifted$ymax, na.rm = TRUE))
  )
  if (shifted$xmin >= shifted$xmax || shifted$ymin >= shifted$ymax) {
    return(NULL)
  }
  shifted
}

coin_df <- load_coin_path(
  coin_path_white_file,
  screen_width,
  screen_height,
  coin_path_screen_w,
  coin_path_screen_h,
  flip_coin_y,
  content_width,
  content_height,
  content_offset_x,
  content_offset_y
)

occ <- load_occluder(occluder_white_file, screen_width, screen_height)
occ <- shift_occluder(occ, white_occluder_dx, white_occluder_dy, screen_width, screen_height)

# ==== Phase config ====
phase_config <- list(
  perception = list(phase_tag = "phase_p", video_col = "video_file_p", trial_col = "trial_p", stage_keep = "0"),
  implicit   = list(phase_tag = "phase_i", video_col = "video_file_i", trial_col = "trial_i", stage_keep = "1"),
  explicit   = list(phase_tag = "phase_e", video_col = "video_file_e", trial_col = "trial_e", stage_keep = "1")
)

parse_video_meta <- function(video_name) {
  base <- tools::file_path_sans_ext(video_name)
  pieces <- strsplit(base, "-")[[1]]
  color <- if (length(pieces)) tolower(pieces[1]) else NA_character_
  stage <- if (length(pieces) >= 2) pieces[2] else "0"
  list(color = color, stage = stage)
}

select_primary_events <- function(df, trial_col) {
  trial_sym <- rlang::sym(trial_col)
  df %>%
    group_by(!!trial_sym, Event) %>%
    summarize(
      samples = n(),
      trial_index = suppressWarnings(first(TRIAL_INDEX)),
      .groups = "drop"
    ) %>%
    group_by(!!trial_sym) %>%
    arrange(trial_index, desc(samples), Event) %>%
    slice(1) %>%
    pull(Event)
}

prepare_gaze <- function(df, trial_col, smooth_gaze, smooth_window) {
  if (!nrow(df)) return(tibble())
  keep_events <- select_primary_events(df, trial_col)
  out <- df %>%
    filter(
      Event %in% keep_events,
      !is.na(Gaze_X),
      !is.na(Gaze_Y),
      is.finite(Time),
      Time >= min_time_ms
    ) %>%
    arrange(Event, Time)
  if (smooth_gaze && smooth_window > 1) {
    out <- out %>%
      group_by(Event) %>%
      mutate(
        Gaze_X = dplyr::coalesce(
          as.numeric(stats::filter(Gaze_X, rep(1 / smooth_window, smooth_window), sides = 2)),
          Gaze_X
        ),
        Gaze_Y = dplyr::coalesce(
          as.numeric(stats::filter(Gaze_Y, rep(1 / smooth_window, smooth_window), sides = 2)),
          Gaze_Y
        )
      ) %>%
      ungroup()
  }
  out
}

compute_avg_path <- function(df) {
  if (!nrow(df)) return(tibble())
  df %>%
    arrange(Event, Time) %>%
    group_by(Event) %>%
    mutate(sample_order = row_number()) %>%
    ungroup() %>%
    group_by(sample_order) %>%
    summarize(
      Gaze_X = mean(Gaze_X, na.rm = TRUE),
      Gaze_Y = mean(Gaze_Y, na.rm = TRUE),
      Time = mean(Time, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    drop_na(Gaze_X, Gaze_Y) %>%
    arrange(sample_order)
}

line_colors <- c(
  perception = "#E76F51",
  explicit = "#2A9D8F",
  implicit = "#457B9D"
)

phase_labels <- c(
  perception = "Perception",
  explicit = "Instructed Imagery",
  implicit = "Passive Viewing"
)

make_plot <- function(participant_id) {
  rds_file <- file.path(rds_dir, paste0("Data_", participant_id, ".rds"))
  if (!file.exists(rds_file)) {
    warning("RDS not found: ", rds_file)
    return(NULL)
  }
  data <- readRDS(rds_file)
  if (!nrow(data)) return(NULL)

  phase_paths <- list()
  for (phase_label in names(phase_config)) {
    cfg <- phase_config[[phase_label]]
    phase_data <- data %>%
      filter(phase == cfg$phase_tag) %>%
      mutate(
        video_name = .data[[cfg$video_col]],
        trial_id = .data[[cfg$trial_col]]
      ) %>%
      filter(!is.na(video_name), video_name != ".", !is.na(trial_id))
    if (!nrow(phase_data)) next

    meta_color <- vapply(phase_data$video_name, function(x) parse_video_meta(x)$color, character(1))
    meta_stage <- vapply(phase_data$video_name, function(x) parse_video_meta(x)$stage, character(1))
    phase_data <- phase_data %>%
      mutate(color = meta_color, stage = meta_stage) %>%
      filter(color == "white", stage == cfg$stage_keep)
    if (!nrow(phase_data)) next

    gaze_df <- prepare_gaze(phase_data, cfg$trial_col, smooth_gaze, smooth_window)
    if (!nrow(gaze_df)) next
    avg_path <- compute_avg_path(gaze_df)
    if (!nrow(avg_path)) next
    avg_path$Participant <- participant_id
    avg_path$Source <- phase_label
    phase_paths[[phase_label]] <- avg_path
  }

  if (!length(phase_paths)) return(NULL)
  plot_df <- bind_rows(phase_paths)
  plot_df$Source <- factor(plot_df$Source, levels = c("perception", "explicit", "implicit"))

  p <- ggplot() +
    {if (!is.null(occ) && nrow(occ)) geom_rect(
      data = occ,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "#BDBDBD",
      alpha = 0.2,
      color = "#666666",
      linewidth = 0.6
    )} +
    geom_path(
      data = plot_df,
      aes(x = Gaze_X, y = Gaze_Y, color = Source),
      linewidth = 1.1
    ) +
    scale_color_manual(values = line_colors, labels = phase_labels) +
    scale_x_continuous(limits = c(0, screen_width), labels = comma) +
    scale_y_reverse(limits = c(screen_height, 0), labels = comma) +
    coord_equal() +
    labs(
      title = paste0("Participant ", participant_id),
      x = "Gaze X (px)",
      y = "Gaze Y (px)",
      color = "Phase"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      text = element_text(family = "Times New Roman")
    )

  if (!is.null(coin_df) && nrow(coin_df)) {
    p <- p +
      geom_path(
        data = coin_df,
        aes(x = screen_x, y = screen_y),
        inherit.aes = FALSE,
        color = "#F2C94C",
        linewidth = 1,
        linetype = "dashed",
        alpha = 0.85
      )
  }
  p
}

plots <- lapply(participants, make_plot)
plots <- plots[!vapply(plots, is.null, logical(1))]
if (length(plots) == 0) stop("No plots generated.")

png(output_path, width = 18, height = 6.6, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, length(plots), heights = unit(c(1, 0.07), "null"))))
for (i in seq_along(plots)) {
  print(plots[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i))
}
grid.text(
  caption_text,
  vp = viewport(layout.pos.row = 2, layout.pos.col = 1:length(plots)),
  gp = gpar(fontfamily = "Times New Roman", fontsize = 11)
)
dev.off()

cat("Saved:", output_path, "\n")
