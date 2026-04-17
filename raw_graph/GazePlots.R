library(dplyr)
library(ggplot2)
library(rlang)
library(scales)
library(tibble)
library(purrr)
library(tidyr)
library(readxl)

# input path
default_data_dir <- "data/processed/test_input"
save_avg_plots <- FALSE
save_color_pairs <- FALSE
merge_in_root <- TRUE
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  data_dir <- normalizePath(args[1], mustWork = TRUE)
} else {
  data_dir <- default_data_dir
}

rds_guess <- list.files(data_dir, pattern = "^Data_.*\\.rds$", full.names = TRUE)
if (length(args) >= 2) {
  rds_file <- normalizePath(args[2], mustWork = TRUE)
} else if (length(rds_guess) == 1) {
  rds_file <- rds_guess
} else if (length(rds_guess) > 1) {
  stop("multiple RDS files：", paste(rds_guess, collapse = ", "))
} else {
  rds_file <- file.path(data_dir, paste0("Data_", basename(data_dir), ".rds"))
}

if (!file.exists(rds_file)) {
  stop("no RDS file：", rds_file)
}

data <- readRDS(rds_file)
if (!nrow(data)) stop("no RDS file")

# ==== Participant Label ====
dir_name <- basename(data_dir)
participant_name <- sub("_.*$", "", dir_name)
if (participant_name %in% c("data", "data_test", "overall", "", basename(dirname(data_dir)))) {
  rds_base <- tools::file_path_sans_ext(basename(rds_file))
  participant_name <- sub("^Data_", "", rds_base)
  participant_name <- sub("_.*$", "", participant_name)
}
participant_name <- toupper(participant_name)
cat("👤 Participant ID detected:", participant_name, "\n")


dir_avg <- if (save_avg_plots) file.path(data_dir, "GazePlots_avg") else NULL
dir_merge <- if (merge_in_root) data_dir else file.path(data_dir, "GazePlots_colorMerged")
dir_pairs <- if (save_color_pairs) file.path(data_dir, "GazePlots_colorPairs") else NULL
if (!is.null(dir_avg)) dir.create(dir_avg, showWarnings = FALSE, recursive = TRUE)
if (!is.null(dir_merge)) dir.create(dir_merge, showWarnings = FALSE, recursive = TRUE)
if (!is.null(dir_pairs)) dir.create(dir_pairs, showWarnings = FALSE, recursive = TRUE)
if (!is.null(dir_avg)) cat("📂 Averaged gaze plots ->", dir_avg, "\n")
cat("📂 Phase-color merged plots ->", dir_merge, "\n")
if (!is.null(dir_pairs)) cat("📂 Color pair plots ->", dir_pairs, "\n")

# Screen size for plotting (pixels)
screen_width <- 4480
screen_height <- 2520

# ==== Coin path (Unity export) ====
coin_path_white_file <- "data/reference/coin_path_white.xlsx"
coin_path_black_file <- "data/reference/coin_path_black.xlsx"
coin_path_screen_w <- screen_width
coin_path_screen_h <- screen_height
occluder_white_file <- "data/reference/white_1_occluder.csv"
occluder_black_file <- "data/reference/black_1_occluder.csv"
white_occluder_dx <- 0
white_occluder_dy <- -200
black_occluder_dx <- -400
black_occluder_dy <- 0

load_coin_path <- function(path, screen_width, screen_height, base_w, base_h) {
  if (!file.exists(path)) return(NULL)
  df <- readxl::read_excel(path)
  if (!all(c("screen_x", "screen_y") %in% names(df))) {
    warning("Coin path no screen_x/screen_y: ", path)
    return(NULL)
  }
  if ("t" %in% names(df)) {
    df <- df %>% arrange(t)
  }
  base_w <- as.numeric(base_w)
  base_h <- as.numeric(base_h)
  if (is.na(base_w) || base_w <= 0) base_w <- screen_width
  if (is.na(base_h) || base_h <= 0) base_h <- screen_height
  sx <- screen_width / base_w
  sy <- screen_height / base_h
  df %>%
    transmute(
      screen_x = as.numeric(screen_x) * sx,
      screen_y = as.numeric(screen_y) * sy
    ) %>%
    drop_na(screen_x, screen_y)
}

coin_paths <- list(
  white = load_coin_path(coin_path_black_file, screen_width, screen_height, coin_path_screen_w, coin_path_screen_h),
  black = load_coin_path(coin_path_white_file, screen_width, screen_height, coin_path_screen_w, coin_path_screen_h)
)

load_occluder <- function(path, screen_width, screen_height) {
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("left", "right", "bottom", "top") %in% names(df))) {
    warning("Occluder is missing left/right/bottom/top: ", path)
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
  # Unity axis opposite
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
  data.frame(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax
  )
}

occluders <- list(
  white = load_occluder(occluder_white_file, screen_width, screen_height),
  black = load_occluder(occluder_black_file, screen_width, screen_height)
)
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
occluders$white <- shift_occluder(occluders$white, white_occluder_dx, white_occluder_dy, screen_width, screen_height)
occluders$black <- shift_occluder(occluders$black, black_occluder_dx, black_occluder_dy, screen_width, screen_height)


phase_config <- list(
  perception = list(phase_tag = "phase_p", video_col = "video_file_p", trial_col = "trial_p"),
  implicit   = list(phase_tag = "phase_i", video_col = "video_file_i", trial_col = "trial_i"),
  explicit   = list(phase_tag = "phase_e", video_col = "video_file_e", trial_col = "trial_e"),
  OE         = list(phase_tag = "phase_oe", video_col = "video_file_oe", trial_col = "trial_oe")
)

required_cols <- unique(unlist(lapply(phase_config, function(cfg) c("phase", cfg$video_col, cfg$trial_col))))
required_cols <- c(required_cols, "Event", "Time", "Gaze_X", "Gaze_Y")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols)) {
  stop("no data：", paste(missing_cols, collapse = ", "))
}

select_primary_events <- function(df, trial_col) {
  trial_sym <- sym(trial_col)
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

prepare_gaze <- function(df, trial_col) {
  if (!nrow(df)) return(tibble())
  keep_events <- select_primary_events(df, trial_col)
  trial_sym <- sym(trial_col)
  df %>%
    filter(Event %in% keep_events, !is.na(Gaze_X), !is.na(Gaze_Y)) %>%
    mutate(trial_factor = factor(.data[[trial_col]])) %>%
    arrange(Event, Time)
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


parse_video_meta <- function(video_name) {
  base <- tools::file_path_sans_ext(video_name)
  pieces <- strsplit(base, "-")[[1]]
  color <- if (length(pieces)) tolower(pieces[1]) else "unknown"
  stage <- if (length(pieces) >= 2) pieces[2] else "0"
  list(base = base, color = color, stage = stage)
}


average_paths <- function(path_list) {
  bind_rows(path_list, .id = "source") %>%
    group_by(sample_order) %>%
    summarize(
      Gaze_X = mean(Gaze_X, na.rm = TRUE),
      Gaze_Y = mean(Gaze_Y, na.rm = TRUE),
      Time = mean(Time, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(sample_order)
}

avg_catalog <- tibble()


for (phase_label in names(phase_config)) {
  cfg <- phase_config[[phase_label]]
  phase_data <- data %>% filter(phase == cfg$phase_tag)
  videos <- phase_data[[cfg$video_col]]
  videos <- unique(videos[!is.na(videos) & videos != "."])
  if (!length(videos)) {
    cat("ℹ️" , phase_label, "skip\n")
    next
  }
  cat("\n======", phase_label, "======\n")
  for (video_name in videos) {
    df <- phase_data %>% filter(.data[[cfg$video_col]] == video_name)
    gaze_df <- prepare_gaze(df, cfg$trial_col)
    if (!nrow(gaze_df)) {
      cat("no gaze data:", phase_label, video_name, "\n")
      next
    }
    n_trials <- dplyr::n_distinct(gaze_df$trial_factor)
    avg_path <- compute_avg_path(gaze_df)
    if (!nrow(avg_path)) {
      cat("no data:", phase_label, video_name, "\n")
      next
    }
    video_meta <- parse_video_meta(video_name)
    color_palette <- c(white = "#F2C94C", black = "#1E1E1E")
    line_color <- color_palette[[video_meta$color]]
    if (is.null(line_color)) line_color <- "#2C7FB8"
    title_text <- paste0(tools::toTitleCase(phase_label), " gaze trajectory (", participant_name, ")")
    subtitle_text <- paste0(video_name, " · mean of ", n_trials, " trials")
    plot_path <- ggplot(avg_path, aes(x = Gaze_X, y = Gaze_Y)) +
      geom_path(color = line_color, linewidth = 1) +
      scale_x_continuous(limits = c(0, screen_width), labels = comma) +
      scale_y_reverse(limits = c(screen_height, 0), labels = comma) +
      coord_equal() +
      labs(
        title = title_text,
        subtitle = subtitle_text,
        x = "Gaze X (px)",
        y = "Gaze Y (px)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "none"
      )
    file_name <- sprintf("gaze2d_%s_%s_%s.png", phase_label, tools::file_path_sans_ext(video_name), participant_name)
    if (save_avg_plots && !is.null(dir_avg)) {
      ggsave(file.path(dir_avg, file_name), plot_path, width = 6, height = 5, dpi = 300)
      cat("saved:", file.path(dir_avg, file_name), "\n")
    }

    avg_catalog <- bind_rows(
      avg_catalog,
      tibble(
        phase = phase_label,
        video = video_name,
        color = video_meta$color,
        stage = video_meta$stage,
        n_trials = n_trials,
        path = list(avg_path)
      )
    )
  }
}

if (!nrow(avg_catalog)) {
  stop("No averaged trajectories are available for merging.")
}

# Combine Perception + Implicit/Explicit
stage_targets <- sort(unique(avg_catalog$stage[avg_catalog$phase %in% c("implicit", "explicit")]))
color_ids <- c("white", "black")

stage_color_components <- list()

for (stage_id in stage_targets) {
  for (color_id in color_ids) {
    components <- avg_catalog %>%
      filter(
        color == color_id,
        (phase == "perception" & stage == "0") |
          (phase %in% c("implicit", "explicit") & stage == stage_id)
      )
    if (!nrow(components)) next
    plot_df <- purrr::map2_dfr(
      components$path,
      seq_len(nrow(components)),
      ~mutate(.x, Source = paste0(tools::toTitleCase(components$phase[.y]), " · ", tools::file_path_sans_ext(components$video[.y])))
    )
    if (!nrow(plot_df)) next
    stage_color_components[[stage_id]][[color_id]] <- list(
      components = components,
      plot_df = plot_df
    )
    plot_stage <- ggplot(plot_df, aes(x = Gaze_X, y = Gaze_Y, color = Source))
    if (stage_id == "1") {
      occ <- occluders[[color_id]]
      if (!is.null(occ) && nrow(occ)) {
        plot_stage <- plot_stage +
          geom_rect(
            data = occ,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#8C8C8C",
            alpha = 0.22,
            color = NA
          ) +
          geom_rect(
            data = occ,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = NA,
            color = "#666666",
            linewidth = 0.8,
            linetype = "dotdash",
            alpha = 0.8
          )
      }
    }
    plot_stage <- plot_stage +
      geom_path(linewidth = 1) +
      scale_color_brewer(palette = "Set2") +
      scale_x_continuous(limits = c(0, screen_width), labels = comma) +
      scale_y_reverse(limits = c(screen_height, 0), labels = comma) +
      coord_equal() +
      labs(
        title = paste0("Stage ", stage_id, " ", color_id, " (", participant_name, ")"),
        subtitle = "Perception + Implicit + Explicit",
        x = "Gaze X (px)",
        y = "Gaze Y (px)",
        color = "Source"
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"), legend.position = "top")
    coin_df <- coin_paths[[color_id]]
    if (!is.null(coin_df) && nrow(coin_df)) {
      coin_color <- if (color_id == "white") "#F2C94C" else "#1E1E1E"
      plot_stage <- plot_stage +
        geom_path(
          data = coin_df,
          aes(x = screen_x, y = screen_y),
          inherit.aes = FALSE,
          color = coin_color,
          linewidth = 1,
          linetype = "dashed",
          alpha = 0.8
        )
    }
    merged_name <- sprintf("gaze_stage%s_%s_%s.png", stage_id, color_id, participant_name)
    ggsave(file.path(dir_merge, merged_name), plot_stage, width = 6, height = 5, dpi = 300)
    cat("Saved merged stage plot:", file.path(dir_merge, merged_name), "\n")
  }
}

# black + white in same graph
for (stage_id in names(stage_color_components)) {
  stage_entry <- stage_color_components[[stage_id]]
  if (!all(c("white", "black") %in% names(stage_entry))) next
  combined_df <- bind_rows(
    stage_entry$black$plot_df %>% mutate(Stimulus = "Black"),
    stage_entry$white$plot_df %>% mutate(Stimulus = "White")
  )
  plot_pair <- ggplot(combined_df, aes(x = Gaze_X, y = Gaze_Y, color = interaction(Stimulus, Source, drop = TRUE))) +
    geom_path(linewidth = 1.1) +
    scale_color_brewer(palette = "Dark2") +
    scale_x_continuous(limits = c(0, screen_width), labels = comma) +
    scale_y_reverse(limits = c(screen_height, 0), labels = comma) +
    coord_equal() +
    labs(
      title = paste0("Stage ", stage_id, " white vs black (", participant_name, ")"),
      x = "Gaze X (px)",
      y = "Gaze Y (px)",
      color = "Stimulus · Source"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  pair_name <- sprintf("gaze_stage%s_black_white_%s.png", stage_id, participant_name)
  if (save_color_pairs && !is.null(dir_pairs)) {
    ggsave(file.path(dir_pairs, pair_name), plot_pair, width = 6, height = 5, dpi = 300)
    cat("Saved black-vs-white plot:", file.path(dir_pairs, pair_name), "\n")
  }
}

cat("\nAveraged and merged gaze graphs finished.\n")
