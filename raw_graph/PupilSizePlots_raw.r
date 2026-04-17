library(ggplot2) 
library(dplyr) 

# input
default_input <- "data/processed/RDS_Data" 
default_output_root <- "outputs/PupilSize_Raw" 

args <- commandArgs(trailingOnly = TRUE) 
if (length(args) >= 1) {
  input_path <- normalizePath(args[0 + 1], mustWork = TRUE) 
} else {
  input_path <- default_input 
}
if (length(args) >= 2) {
  output_root <- normalizePath(args[0 + 2], mustWork = FALSE) 
} else {
  output_root <- default_output_root 
}
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE) 

if (dir.exists(input_path)) {
  rds_files <- list.files(input_path, pattern = "^Data_.*\\.rds$", full.names = TRUE, ignore.case = TRUE) 
  rds_files <- sort(rds_files) 
} else if (file.exists(input_path)) {
  rds_files <- input_path 
} else {
  stop("no file：", input_path) 
}

if (!length(rds_files)) {
  stop("no file：", input_path) 
}

# Parse each video filename into stimulus and suffix labels.
parse_video_info <- function(video_name) {
  base <- tools::file_path_sans_ext(basename(video_name)) 
  stimulus <- if (grepl("white", base, ignore.case = TRUE)) {
    "white"
  } else if (grepl("black", base, ignore.case = TRUE)) {
    "black"
  } else {
    NA_character_
  }
  suffix <- sub(".*-(\\d+)$", "\\1", base) 
  if (identical(suffix, base)) {
    suffix <- NA_character_ 
  }
  list(stimulus = stimulus, suffix = suffix) 
}

get_shadow_window <- function(suffix_label) {
  if (suffix_label == "1") return(c(8000, 12000)) 
  if (suffix_label == "2") return(c(15000, 17000)) 
  if (suffix_label == "3") return(c(20000, 24000)) 
  NULL 
}

# Clean and align raw samples for a single phase-video combination.
prepare_video_trials <- function(df, trial_col, phase_label, video_name) {
  if (!nrow(df)) {
    cat("skip: no data ->", phase_label, video_name, "\n") 
    return(NULL) 
  }

  df <- df %>% 
    filter(!is.na(.data[[trial_col]]), !is.na(Event))

  trial_event_stats <- df %>% 
    group_by(.data[[trial_col]], Event) %>%
    summarize(
      samples = n(),
      trial_index = suppressWarnings(first(TRIAL_INDEX)),
      .groups = "drop"
    )

  primary_events <- trial_event_stats %>% 
    group_by(.data[[trial_col]]) %>%
    arrange(trial_index, desc(samples), Event) %>%
    slice(1) %>%
    ungroup()

  if (nrow(trial_event_stats) > nrow(primary_events)) {
    cat("Info:", phase_label, video_name, ": multiple events found in the same trial; keeping the one with the most samples.\n") 
  }

  df <- df %>% 
    semi_join(primary_events, by = "Event")

  df <- df %>% 
    arrange(Event, Time) %>%
    group_by(Event) %>%
    mutate(time_rel = Time - min(Time, na.rm = TRUE)) %>%
    ungroup()

  if (!nrow(df)) {
    cat("Skip: zero valid rows ->", phase_label, video_name, "\n") 
    return() 
  }

  df <- df %>% 
    mutate(
      trial_factor = factor(.data[[trial_col]]),
      legend_label = as.character(trial_factor)
    )

  if (all(is.na(df$trial_factor))) {
    cat("Skip: trial grouping is empty ->", phase_label, video_name, "\n") 
    return(NULL) 
  }

  df 
}

# Build and save pupil traces for one phase and suffix across all matching videos.
save_pupil_plot <- function(df, trial_col, video_col, phase_label, participant_name, output_dir, suffix_label) {
  if (!nrow(df)) {
    cat("Skip: no data ->", phase_label, suffix_label, "\n") 
    return()
  }

  videos <- unique(df[[video_col]]) 
  videos <- videos[!is.na(videos) & videos != "."] 
  if (!length(videos)) {
    cat("Skip: no valid videos ->", phase_label, suffix_label, "\n") 
    return() 
  }

  combined <- list() 
  for (video_name in videos) {
    df_video <- df %>% 
      filter(.data[[video_col]] == video_name)
    df_video <- prepare_video_trials(df_video, trial_col, phase_label, video_name) 
    if (is.null(df_video) || !nrow(df_video)) next 

    info <- parse_video_info(video_name) 
    if (is.na(info$stimulus)) {
      cat("Skip: could not identify video stimulus ->", phase_label, video_name, "\n") 
      next 
    }
    df_video <- df_video %>% 
      mutate(
        stimulus = info$stimulus,
        group_id = interaction(video_name, trial_factor, drop = TRUE)
      )
    combined[[length(combined) + 1]] <- df_video 
  }

  if (!length(combined)) {
    cat("Skip: no plottable data ->", phase_label, suffix_label, "\n") 
    return() 
  }

  df_all <- bind_rows(combined) 
  df_all <- df_all %>% 
    mutate(stimulus = factor(stimulus, levels = c("white", "black")))


  x_limit <- 30000 
  x_breaks <- seq(0, x_limit, by = 5000) 

  
  y_limits <- c(-4, 4) 
  y_breaks <- seq(-4, 4, by = 2) 
  shadow_window <- get_shadow_window(suffix_label) 

  plot_title <- paste0( 
    tools::toTitleCase(phase_label),
    " Phase_", suffix_label, "_", participant_name
  )


  p <- ggplot( 
    df_all,
    aes(
      x = time_rel,
      y = Pupil,
      color = stimulus,
      group = group_id
    )
  )
  if (!is.null(shadow_window)) {
    p <- p + annotate( 
      "rect",
      xmin = shadow_window[1],
      xmax = shadow_window[2],
      ymin = -Inf,
      ymax = Inf,
      fill = "grey80",
      alpha = 0.25
    )
  }
  p <- p + 
    geom_vline(xintercept = c(7000, 25000), color = "black", linewidth = 0.3) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    scale_x_continuous(
      limits = c(0, x_limit),
      breaks = x_breaks,
      labels = function(x) format(x, big.mark = ",", trim = TRUE)
    ) +
    scale_y_continuous(limits = range(y_limits), breaks = y_breaks) +
    scale_color_manual(values = c(white = "#FFD60A", black = "#0B3D91")) +
    labs(
      title = plot_title,
      x = "Time (ms)",
      y = "Pupil size (mm)",
      color = "Stimulus"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )

  file_name <- paste0(phase_label, "_bw-", suffix_label, "_", participant_name, ".png") 
  file_path <- file.path(output_dir, file_name) 
  ggsave(file_path, plot = p, width = 8, height = 4.5, dpi = 300) 
  cat("save：", file_path, "\n") 
}


phase_config <- list( 
  perception = list(phase_tag = "phase_p", video_col = "video_file_p", trial_col = "TRIAL_INDEX"),
  implicit   = list(phase_tag = "phase_i", video_col = "video_file_i", trial_col = "TRIAL_INDEX"),
  explicit   = list(phase_tag = "phase_e", video_col = "video_file_e", trial_col = "TRIAL_INDEX"),
  OE         = list(phase_tag = "phase_oe", video_col = "video_file_oe", trial_col = "TRIAL_INDEX")
)

# Process one participant RDS and export phase-by-suffix pupil plots.
process_one_rds <- function(rds_file, output_root) {
  cat("\nProcessing RDS:", rds_file, "\n") 
  data <- readRDS(rds_file) 

  if (!nrow(data)) {
    cat("skip：", rds_file, "\n") 
    return(invisible(FALSE)) 
  }

  required_cols <- c("phase", "Event") 
  missing_cols <- setdiff(required_cols, names(data)) 
  if (length(missing_cols)) {
    cat("Warning: missing required columns, skipping:", rds_file, "-", paste(missing_cols, collapse = ", "), "\n")
    return(invisible(FALSE)) 
  }

  participant_name <- basename(rds_file) 
  participant_name <- sub("^Data_", "", participant_name) 
  participant_name <- sub("\\.[Rr][Dd][Ss]$", "", participant_name) 

  output_dir <- file.path(output_root, participant_name) 
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE) 
  cat("Participant ID detected:", participant_name, "\n") 
  cat("save：", output_dir, "\n") 


  for (phase_label in names(phase_config)) {
    cfg <- phase_config[[phase_label]] 
    required_cols <- c("phase", cfg$video_col, cfg$trial_col) 
    missing_cols <- setdiff(required_cols, names(data)) 
    if (length(missing_cols)) {
      cat("skip：", phase_label, "-", paste(missing_cols, collapse = ", "), "\n") 
      next 
    }

    phase_data <- data %>% 
      filter(phase == cfg$phase_tag)

    videos <- phase_data[[cfg$video_col]] 
    videos <- unique(videos[!is.na(videos) & videos != "."]) 

    if (!length(videos)) {
      cat("no usable video：", phase_label, "\n") 
      next 
    }

    video_info <- lapply(videos, parse_video_info) 
    video_table <- data.frame( 
      video = videos,
      stimulus = vapply(video_info, function(x) x$stimulus, character(1)),
      suffix = vapply(video_info, function(x) x$suffix, character(1)),
      stringsAsFactors = FALSE
    )
    video_table <- video_table[!is.na(video_table$stimulus) & !is.na(video_table$suffix), ] 
    if (!nrow(video_table)) {
      cat("no usable video：", phase_label, "\n") 
      next 
    }

    cat("\n====== Processing", phase_label, "phase ======\n") 

    for (suffix_label in unique(video_table$suffix)) {
      group_videos <- video_table$video[video_table$suffix == suffix_label] 
      df <- phase_data %>% 
        filter(.data[[cfg$video_col]] %in% group_videos)

      save_pupil_plot(df, cfg$trial_col, cfg$video_col, phase_label, participant_name, output_dir, suffix_label) 
    }
  }

  cat("\n saved：", output_dir, "\n") 
  invisible(TRUE) 
}

for (rds_file in rds_files) {
  process_one_rds(rds_file, output_root) 
}
