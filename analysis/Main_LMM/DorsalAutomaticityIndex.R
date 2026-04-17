library(dplyr) 
library(readxl) 
library(tidyr) 


default_input <- "data/processed/RDS_Data" 
default_output_dir <- "outputs/PupilSize_Raw" 

coin_path_white_file <- "data/reference/White_Trajectory.xlsx" 
coin_path_black_file <- "data/reference/Black_Trajectory.xlsx" 
swap_by_phase <- list(perception = FALSE, implicit = FALSE, explicit = FALSE)


screen_width <- 4480 
screen_height <- 2520 

content_width <- 3840
content_height <- 2160
content_offset_x <- 0
content_offset_y <- 0


time_start_ms <- 8000 
time_end_ms <- 12000 
time_shift_ms <- 2000 

args <- commandArgs(trailingOnly = TRUE) 
if (length(args) >= 1) {
  input_path <- normalizePath(args[1], mustWork = TRUE) 
} else {
  input_path <- default_input 
}
if (length(args) >= 2) {
  output_dir <- normalizePath(args[2], mustWork = FALSE) 
} else {
  output_dir <- default_output_dir 
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE) 

trial_output <- file.path(output_dir, "dorsal_trial_errors.csv") 

if (dir.exists(input_path)) {
  rds_files <- list.files(input_path, pattern = "^Data_.*\\.rds$", full.names = TRUE, ignore.case = TRUE) 
  rds_files <- sort(rds_files) 
} else if (file.exists(input_path)) {
  rds_files <- input_path 
} else {
  stop("no RDS file：", input_path) 
}

if (!length(rds_files)) {
  stop("no RDS file：", input_path) 
}


load_coin_path <- function(path, screen_width, screen_height,
                           content_width, content_height, offset_x, offset_y) {
  if (!file.exists(path)) return(NULL) 
  df <- readxl::read_excel(path)
  has_viewport <- all(c("viewport_x", "viewport_y") %in% names(df))
  has_screen <- all(c("screen_x", "screen_y") %in% names(df))
  if (!has_viewport && !has_screen) {
    stop("Coin path is missing viewport_x/viewport_y or screen_x/screen_y: ", path)
  }
  if ("t" %in% names(df)) {
    df <- df %>% arrange(t)
  }
  if (has_viewport) {
    out <- df %>%
      transmute(
        t_ms = as.numeric(t) * 1000, 
        screen_x = as.numeric(viewport_x) * content_width + offset_x,
        screen_y = (1 - as.numeric(viewport_y)) * content_height + offset_y
      ) %>%
      drop_na(t_ms, screen_x, screen_y)
  } else {
    out <- df %>%
      transmute(
        t_ms = as.numeric(t) * 1000,
        screen_x = as.numeric(screen_x),
        screen_y = as.numeric(screen_y)
      ) %>%
      drop_na(t_ms, screen_x, screen_y)
  }
  out
}

prepare_target_path <- function(path_df) {
  if (is.null(path_df) || !nrow(path_df)) return(NULL) 
  x <- as.numeric(path_df$screen_x)
  y <- as.numeric(path_df$screen_y)
  t_ms <- as.numeric(path_df$t_ms)
  dx <- diff(x)
  dy <- diff(y)
  seg_len <- sqrt(dx * dx + dy * dy)
  list(
    t_ms = t_ms,
    x = x,
    y = y,
    dx = dx,
    dy = dy,
    seg_len = seg_len,
    seg_len2 = dx * dx + dy * dy,
    cum_len = c(0, cumsum(seg_len))
  )
}

interp_target <- function(t_ms, path) {
  list(
    x = approx(path$t_ms, path$x, xout = t_ms, ties = mean, rule = 2)$y,
    y = approx(path$t_ms, path$y, xout = t_ms, ties = mean, rule = 2)$y,
    s = approx(path$t_ms, path$cum_len, xout = t_ms, ties = mean, rule = 2)$y
  )
}

project_onto_path <- function(xg, yg, path) {
  if (is.na(xg) || is.na(yg)) return(c(s_proj = NA_real_, dist = NA_real_))
  x1 <- path$x[-length(path$x)]
  y1 <- path$y[-length(path$y)]
  dx <- path$dx
  dy <- path$dy
  seg_len2 <- path$seg_len2
  valid <- seg_len2 > 0
  if (!any(valid)) return(c(s_proj = NA_real_, dist = NA_real_))
  t <- ((xg - x1) * dx + (yg - y1) * dy) / seg_len2
  t <- pmax(0, pmin(1, t))
  proj_x <- x1 + t * dx
  proj_y <- y1 + t * dy
  dist2 <- (xg - proj_x) * (xg - proj_x) + (yg - proj_y) * (yg - proj_y)
  dist2[!valid] <- Inf
  idx <- which.min(dist2)
  c(
    s_proj = path$cum_len[idx] + path$seg_len[idx] * t[idx],
    dist = sqrt(dist2[idx])
  )
}

compute_trial_metrics <- function(trial_df, target_path, time_shift_ms = 0) {
  if (!nrow(trial_df)) return(list(rmse = NA_real_, n = 0)) 

  target_interp <- interp_target(trial_df$Time - time_shift_ms, target_path) 
  valid <- !is.na(trial_df$Gaze_X) &
    !is.na(trial_df$Gaze_Y) &
    !is.na(target_interp$x) &
    !is.na(target_interp$y)
  if (!any(valid)) return(list(rmse = NA_real_, n = 0))
  gaze_x <- trial_df$Gaze_X[valid]
  gaze_y <- trial_df$Gaze_Y[valid]
  dx <- gaze_x - target_interp$x[valid]
  dy <- gaze_y - target_interp$y[valid]
  dist_target <- sqrt(dx^2 + dy^2) 
  rmse <- if (!length(dist_target)) NA_real_ else sqrt(mean(dist_target^2, na.rm = TRUE))
  list(
    rmse = rmse,
    n = sum(valid)
  )
}

extract_stimulus <- function(video_name) {
  if (is.na(video_name)) return(NA_character_)
  if (grepl("white", video_name, ignore.case = TRUE)) return("white")
  if (grepl("black", video_name, ignore.case = TRUE)) return("black")
  NA_character_
}

extract_suffix <- function(video_name) {
  if (is.na(video_name)) return(NA_character_)
  base <- tools::file_path_sans_ext(basename(video_name))
  suffix <- sub(".*-(\\d+)$", "\\1", base)
  if (identical(suffix, base)) return(NA_character_)
  suffix
}

pad_id <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.integer(x_chr))
  ifelse(is.na(x_num), x_chr, sprintf("%08d", x_num))
}


coin_white <- load_coin_path(
  coin_path_white_file,
  screen_width,
  screen_height,
  content_width,
  content_height,
  content_offset_x,
  content_offset_y
)
coin_black <- load_coin_path(
  coin_path_black_file,
  screen_width,
  screen_height,
  content_width,
  content_height,
  content_offset_x,
  content_offset_y
)

build_targets <- function(swapped) {
  if (swapped) {
    list(
      white = prepare_target_path(coin_black),
      black = prepare_target_path(coin_white)
    )
  } else {
    list(
      white = prepare_target_path(coin_white),
      black = prepare_target_path(coin_black)
    )
  }
}

targets_default <- build_targets(FALSE)
targets_swapped <- build_targets(TRUE)

if (is.null(targets_default$white) || is.null(targets_default$black)) {
  stop("no coin path")
}


phase_config <- list(
  perception = list(phase_tag = "phase_p", video_col = "video_file_p", trial_col = "TRIAL_INDEX"),
  implicit = list(phase_tag = "phase_i", video_col = "video_file_i", trial_col = "TRIAL_INDEX"),
  explicit = list(phase_tag = "phase_e", video_col = "video_file_e", trial_col = "TRIAL_INDEX")
)

trial_results <- list() 

for (rds_file in rds_files) {
  cat("Processing RDS:", rds_file, "\n") 
  data <- readRDS(rds_file) 
  if (!nrow(data)) {
    cat("no RDS file：", rds_file, "\n") 
    next
  }
  
  required_cols <- c(
    "phase",
    "Time",
    "Gaze_X",
    "Gaze_Y",
    "video_file_p",
    "video_file_i",
    "video_file_e",
    "TRIAL_INDEX"
  ) 
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols)) {
    cat("skip：", rds_file, "-", paste(missing_cols, collapse = ", "), "\n")
    next
  }
  
  participant_name <- basename(rds_file) 
  participant_name <- sub("^Data_", "", participant_name) 
  participant_name <- sub("\\.[Rr][Dd][Ss]$", "", participant_name) 
  participant_name <- pad_id(participant_name)
  
  for (phase_label in names(phase_config)) {
    cfg <- phase_config[[phase_label]] 
    phase_data <- data %>%
      filter(
        phase == cfg$phase_tag,
        Time >= time_start_ms,
        Time <= time_end_ms
      ) %>%
      mutate(
        video_name = .data[[cfg$video_col]],
        trial_id = .data[[cfg$trial_col]]
      ) %>%
      filter(!is.na(video_name), video_name != ".", !is.na(trial_id))
    
    if (!nrow(phase_data)) next
    
    phase_data <- phase_data %>%
      mutate(
        stimulus = vapply(video_name, extract_stimulus, character(1)),
        suffix = vapply(video_name, extract_suffix, character(1))
      )

    if (phase_label == "perception") {
      phase_data <- phase_data %>%
        filter(stimulus %in% c("white", "black"))
    } else {
      phase_data <- phase_data %>%
        filter(suffix == "1", stimulus %in% c("white", "black"))
    }
    
    if (!nrow(phase_data)) next
    
    trial_groups <- phase_data %>%
      group_by(trial_id, video_name, stimulus) %>%
      group_split()
    
    for (trial_df in trial_groups) {
      stim <- unique(trial_df$stimulus)
      if (length(stim) != 1 || is.na(stim)) next
      phase_swap <- isTRUE(swap_by_phase[[phase_label]])
      targets <- if (phase_swap) targets_swapped else targets_default
      target_path <- targets[[stim]]
      if (is.null(target_path)) next
      metrics <- compute_trial_metrics(trial_df, target_path, time_shift_ms = 0)
      metrics_shift <- compute_trial_metrics(trial_df, target_path, time_shift_ms = time_shift_ms)
      trial_results[[length(trial_results) + 1]] <- data.frame(
        participant = participant_name,
        phase = phase_label,
        stimulus = stim,
        trial_id = unique(trial_df$trial_id),
        video_name = unique(trial_df$video_name),
        traj_rmse = metrics$rmse,
        traj_rmse_shift_2000ms = metrics_shift$rmse,
        n_samples = metrics$n,
        n_samples_shift = metrics_shift$n,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (!length(trial_results)) {
  stop("no result")
}

trial_df <- bind_rows(trial_results)

write.csv(trial_df, trial_output, row.names = FALSE)
cat("saved:", trial_output, "\n")
