library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)

# input
default_input <- "data/processed/RDS_Data"
default_output_dir <- "outputs/Gaze_mean_ci"

smooth_gaze <- TRUE
smooth_window <- 5
min_time_ms <- 1000

# Screen size
screen_width <- 4480
screen_height <- 2520

# Coin paths (no swap, flip Y)
coin_path_white_file <- "data/reference/White_Trajectory.xlsx"
coin_path_black_file <- "data/reference/Black_Trajectory.xlsx"
flip_coin_y <- TRUE
content_width <- 3840
content_height <- 2160
content_offset_x <- 0
content_offset_y <- 0

# Occluder (stage 1)
occluder_white_file <- "data/reference/white_1_occluder.csv"
occluder_black_file <- "data/reference/black_1_occluder.csv"
white_occluder_dx <- -150
white_occluder_dy <- -200
black_occluder_dx <- -550
black_occluder_dy <- 0
black_occluder_ymin <- 1000
black_occluder_ymax <- 2250

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else default_input
output_dir <- if (length(args) >= 2) normalizePath(args[2], mustWork = FALSE) else default_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}

phase_config <- list(
  perception = list(phase_tag = "phase_p", video_col = "video_file_p", trial_col = "trial_p", stage_keep = "0"),
  implicit   = list(phase_tag = "phase_i", video_col = "video_file_i", trial_col = "trial_i", stage_keep = "1"),
  explicit   = list(phase_tag = "phase_e", video_col = "video_file_e", trial_col = "trial_e", stage_keep = "1")
)

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
    filter(Event %in% keep_events, !is.na(Gaze_X), !is.na(Gaze_Y), Time >= min_time_ms) %>%
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

compute_mean_ci_path <- function(phase_list, phase_label) {
  if (!length(phase_list)) return(NULL)
  phase_df <- bind_rows(phase_list)
  if (!nrow(phase_df)) return(NULL)
  mean_path <- phase_df %>%
    group_by(sample_order) %>%
    summarize(
      mean_x = mean(Gaze_X, na.rm = TRUE),
      mean_y = mean(Gaze_Y, na.rm = TRUE),
      n = n_distinct(Participant),
      .groups = "drop"
    ) %>%
    arrange(sample_order)
  phase_df <- phase_df %>%
    left_join(mean_path, by = "sample_order") %>%
    mutate(dist = sqrt((Gaze_X - mean_x)^2 + (Gaze_Y - mean_y)^2))
  stats_df <- phase_df %>%
    group_by(sample_order) %>%
    summarize(
      mean_x = first(mean_x),
      mean_y = first(mean_y),
      n = n_distinct(Participant),
      sd_dist = sd(dist, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(sample_order) %>%
    mutate(
      se_dist = ifelse(n > 1, sd_dist / sqrt(n), NA_real_),
      t_crit = ifelse(n > 1, stats::qt(0.975, df = n - 1), NA_real_),
      ci_radius = se_dist * t_crit
    ) %>%
    mutate(
      dx = dplyr::lead(mean_x, default = dplyr::last(mean_x)) -
        dplyr::lag(mean_x, default = dplyr::first(mean_x)),
      dy = dplyr::lead(mean_y, default = dplyr::last(mean_y)) -
        dplyr::lag(mean_y, default = dplyr::first(mean_y)),
      norm_len = sqrt(dx^2 + dy^2),
      nx = ifelse(norm_len > 0, -dy / norm_len, 0),
      ny = ifelse(norm_len > 0, dx / norm_len, 0),
      upper_x = mean_x + nx * ci_radius,
      upper_y = mean_y + ny * ci_radius,
      lower_x = mean_x - nx * ci_radius,
      lower_y = mean_y - ny * ci_radius
    )
  mean_df <- stats_df %>%
    transmute(sample_order, Gaze_X = mean_x, Gaze_Y = mean_y, Source = phase_label)
  ribbon_df <- stats_df %>%
    filter(is.finite(ci_radius)) %>%
    transmute(sample_order, upper_x, upper_y, lower_x, lower_y, Source = phase_label)
  list(mean_df = mean_df, ribbon_df = ribbon_df)
}

load_coin_path <- function(path, screen_width, screen_height, flip_y,
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
    out <- df %>%
      transmute(
        screen_x = as.numeric(screen_x),
        screen_y = as.numeric(screen_y)
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

coin_paths <- list(
  white = load_coin_path(
    coin_path_white_file, screen_width, screen_height, flip_coin_y,
    content_width, content_height, content_offset_x, content_offset_y
  ),
  black = load_coin_path(
    coin_path_black_file, screen_width, screen_height, flip_coin_y,
    content_width, content_height, content_offset_x, content_offset_y
  )
)

occluders <- list(
  white = load_occluder(occluder_white_file, screen_width, screen_height),
  black = load_occluder(occluder_black_file, screen_width, screen_height)
)
occluders$white <- shift_occluder(occluders$white, white_occluder_dx, white_occluder_dy, screen_width, screen_height)
occluders$black <- shift_occluder(occluders$black, black_occluder_dx, black_occluder_dy, screen_width, screen_height)
if (!is.null(occluders$black) && nrow(occluders$black)) {
  occluders$black$ymin <- max(0, black_occluder_ymin)
  occluders$black$ymax <- min(screen_height, black_occluder_ymax)
}

if (dir.exists(input_path)) {
  rds_files <- list.files(input_path, pattern = "^Data_.*[.]rds$", full.names = TRUE, ignore.case = TRUE)
} else if (file.exists(input_path)) {
  rds_files <- input_path
} else {
  stop("Input not found: ", input_path)
}
if (!length(rds_files)) stop("No RDS files found: ", input_path)

line_colors <- c(explicit = "#2A9D8F", implicit = "#457B9D", perception = "#E76F51")
phase_labels <- c(
  perception = "Perception",
  explicit = "Instructed Imagery",
  implicit = "Passive Viewing"
)

agg_paths <- list(white = list(), black = list())

for (rds_file in rds_files) {
  df <- readRDS(rds_file)
  if (!nrow(df)) next

  participant_id <- NA_character_
  if ("Subject" %in% names(df)) {
    participant_id <- pad_id(na.omit(as.character(df$Subject))[1])
  }
  if (is.na(participant_id) || !nzchar(participant_id)) {
    participant_id <- pad_id(sub("^Data_", "", tools::file_path_sans_ext(basename(rds_file))))
  }

  for (phase_label in names(phase_config)) {
    cfg <- phase_config[[phase_label]]
    phase_df <- df %>%
      filter(.data$phase == cfg$phase_tag) %>%
      mutate(
        trial_id = .data[[cfg$trial_col]],
        video_file = .data[[cfg$video_col]]
      ) %>%
      filter(!is.na(video_file))

    if (!nrow(phase_df)) next

    phase_df <- phase_df %>%
      mutate(video_base = tools::file_path_sans_ext(video_file)) %>%
      separate(video_base, into = c("color", "stage"), sep = "-", fill = "right", extra = "drop") %>%
      mutate(
        color = tolower(color),
        stage = ifelse(is.na(stage), "0", stage)
      )

    for (color_id in c("white", "black")) {
      phase_data <- phase_df %>%
        filter(color == color_id, stage == cfg$stage_keep)
      if (!nrow(phase_data)) next
      gaze_df <- prepare_gaze(phase_data, cfg$trial_col, smooth_gaze, smooth_window)
      if (!nrow(gaze_df)) next
      avg_path <- compute_avg_path(gaze_df)
      if (!nrow(avg_path)) next
      avg_path <- avg_path %>%
        mutate(
          Participant = participant_id,
          phase = phase_label
        )
      if (is.null(agg_paths[[color_id]][[phase_label]])) {
        agg_paths[[color_id]][[phase_label]] <- list(avg_path)
      } else {
        agg_paths[[color_id]][[phase_label]][[length(agg_paths[[color_id]][[phase_label]]) + 1]] <- avg_path
      }
    }
  }
}

plot_mean_ci <- function(color_id, phases, include_occluder, out_name, title_prefix) {
  phase_means <- list()
  phase_ribbons <- list()
  for (phase_label in phases) {
    phase_list <- agg_paths[[color_id]][[phase_label]]
    if (!length(phase_list)) next
    ci_result <- compute_mean_ci_path(phase_list, phase_label)
    if (is.null(ci_result)) next
    phase_means[[phase_label]] <- ci_result$mean_df
    if (nrow(ci_result$ribbon_df)) {
      phase_ribbons[[phase_label]] <- ci_result$ribbon_df
    }
  }
  if (!length(phase_means)) return()
  mean_df <- bind_rows(phase_means)
  mean_df$Source <- factor(mean_df$Source, levels = phases)

  ribbon_poly <- list()
  if (length(phase_ribbons)) {
    for (phase_label in names(phase_ribbons)) {
      rb <- phase_ribbons[[phase_label]]
      if (!nrow(rb)) next
      poly <- bind_rows(
        rb %>%
          arrange(sample_order) %>%
          transmute(Gaze_X = upper_x, Gaze_Y = upper_y, Source = phase_label),
        rb %>%
          arrange(desc(sample_order)) %>%
          transmute(Gaze_X = lower_x, Gaze_Y = lower_y, Source = phase_label)
      )
      ribbon_poly[[phase_label]] <- poly
    }
  }
  ribbon_df <- if (length(ribbon_poly)) bind_rows(ribbon_poly) else tibble()

  coin_df <- coin_paths[[color_id]]
  occ <- occluders[[color_id]]
  p <- ggplot() +
    {if (!is.null(occ) && include_occluder) geom_rect(
      data = occ,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "#BDBDBD",
      alpha = 0.2,
      color = "#666666",
      linewidth = 0.6
    )} +
    {if (nrow(ribbon_df)) geom_polygon(
      data = ribbon_df,
      aes(x = Gaze_X, y = Gaze_Y, fill = Source),
      alpha = 0.18,
      color = NA
    )} +
    geom_path(
      data = mean_df,
      aes(x = Gaze_X, y = Gaze_Y, color = Source),
      linewidth = 1.2
    ) +
    scale_color_manual(values = line_colors, labels = phase_labels) +
    scale_fill_manual(values = line_colors, guide = "none") +
    scale_x_continuous(limits = c(0, screen_width)) +
    scale_y_reverse(limits = c(screen_height, 0)) +
    coord_equal() +
    labs(
      x = "Gaze X (px)",
      y = "Gaze Y (px)",
      color = "Phase"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      text = element_text(family = "Times New Roman")
    )

  if (!is.null(coin_df) && nrow(coin_df)) {
    coin_color <- if (color_id == "white") "#F2C94C" else "#1E1E1E"
    p <- p +
      geom_path(
        data = coin_df,
        aes(x = screen_x, y = screen_y),
        inherit.aes = FALSE,
        color = coin_color,
        linewidth = 1,
        linetype = "dashed",
        alpha = 0.85
      )
  }

  ggsave(file.path(output_dir, out_name), p, width = 6.5, height = 5, dpi = 300)
  cat("Saved:", file.path(output_dir, out_name), "\n")
}

for (color_id in c("white", "black")) {
  plot_mean_ci(
    color_id = color_id,
    phases = c("explicit", "implicit"),
    include_occluder = TRUE,
    out_name = sprintf("gaze_mean_ci_explicit_implicit_%s.png", color_id),
    title_prefix = "Explicit + Implicit mean gaze"
  )
  plot_mean_ci(
    color_id = color_id,
    phases = c("explicit"),
    include_occluder = TRUE,
    out_name = sprintf("gaze_mean_ci_explicit_%s.png", color_id),
    title_prefix = "Explicit mean gaze"
  )
  plot_mean_ci(
    color_id = color_id,
    phases = c("perception"),
    include_occluder = FALSE,
    out_name = sprintf("gaze_mean_ci_perception_%s.png", color_id),
    title_prefix = "Perception mean gaze"
  )
}
