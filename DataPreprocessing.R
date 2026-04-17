rm(list=ls(all=TRUE))  

# install.packages("robustbase")
# install.packages("VWPre")
# install.packages("~/Downloads/PupilPre", repos = NULL, type = "source")

library(data.table)
library(PupilPre)
library(ggplot2)
library(dplyr)

# ==== Phase mapping by TRIAL_INDEX ====
# Fixed phase order: perception (1-4) -> OE (5-8) -> implicit (9-18) -> explicit (19-28)
phase_ranges <- list(
  perception = c(1, 4),
  OE = c(5, 8),
  implicit = c(9, 18),
  explicit = c(19, 28)
)

# Map each trial index to the fixed experiment phase labels.
map_phase_by_trial_index <- function(trial_index) {
  phase <- rep(NA_character_, length(trial_index))
  phase[trial_index >= phase_ranges$perception[1] & trial_index <= phase_ranges$perception[2]] <- "phase_p"
  phase[trial_index >= phase_ranges$OE[1] & trial_index <= phase_ranges$OE[2]] <- "phase_oe"
  phase[trial_index >= phase_ranges$implicit[1] & trial_index <= phase_ranges$implicit[2]] <- "phase_i"
  phase[trial_index >= phase_ranges$explicit[1] & trial_index <= phase_ranges$explicit[2]] <- "phase_e"
  phase
}

# input path
default_input <- "data/raw/Raw_Data_Final"
default_output_dir <- "data/processed/RDS_Data"

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

if (dir.exists(input_path)) {
  files <- list.files(input_path, pattern = "\\.txt$", full.names = TRUE, ignore.case = TRUE)
  files <- sort(files)
} else {
  files <- input_path
}

if (length(files) == 0) {
  stop("No .txt files found in: ", input_path)
}

# Run the full pupil preprocessing pipeline for one raw text file.
process_one <- function(file, output_dir) {
  cat("Processing file:", file, "\n")

  base_name <- basename(file)
  base_name <- sub("^SampleReport_", "", base_name)
  base_name <- sub("\\.[Tt][Xx][Tt]$", "", base_name)
  output_file <- file.path(output_dir, paste0("Data_", base_name, ".rds"))
  if (file.exists(output_file)) {
    cat("Skip existing:", output_file, "\n")
    return(output_file)
  }

  # read data
  rawPupil <- fread(file)

  # data preprocessing
  if (!("TRIAL_INDEX" %in% names(rawPupil))) {
    stop("Missing TRIAL_INDEX column in: ", file)
  }

  rawPupil <- ppl_prep_data(data = rawPupil, Subject = "Session_Name_")
  rawPupil <- rawPupil %>%
    mutate(
      phase = map_phase_by_trial_index(TRIAL_INDEX),
      Event = factor(paste0("T", TRIAL_INDEX))
    )
  if (any(is.na(rawPupil$phase))) {
    cat("Warning: some TRIAL_INDEX values could not be mapped to a phase and were dropped.\n")
    rawPupil <- rawPupil %>% filter(!is.na(phase))
  }


  # remove dots outside of screen range
  rawPupil <- ppl_rm_extra_DVcols(rawPupil)
  rawPupil <- create_time_series(data = rawPupil, Adjust = 0)
  rawPupil <- ppl_select_recorded_eye(data = rawPupil, Recording = "R", WhenLandR = "Right")
  rawPupil <- recode_off_screen(data = rawPupil, ScreenSize = c(4480, 2520))

  # remove blink
  cleanedData <- clean_blink(
    rawPupil,
    BlinkPadding = c(100, 100),
    Delta = 5,
    MaxValueRun = 5,
    NAsAroundRun = c(2, 2),
    LogFile = paste0(tempdir(), "/BlinkCleanupLog.rds")
  )

  rm(rawPupil); gc()

  # remove artifact
  cleanedData <- tryCatch(
    clean_artifact(
      cleanedData,
      MADWindow = 20,
      MADConstant = 4,
      MADPadding = c(100, 100),
      MahaConstant = 1,
      Method = "Robust",
      XandY = TRUE,
      Second = TRUE,
      MaxValueRun = 5,
      NAsAroundRun = c(2, 2),
      LogFile = paste0(tempdir(), "/ArtifactCleanupLog.rds")
    ),
    error = function(e1) {
      cat("WARN: clean_artifact Robust failed, retrying with Mahal. Error:", conditionMessage(e1), "\n")
      tryCatch(
        clean_artifact(
          cleanedData,
          MADWindow = 20,
          MADConstant = 4,
          MADPadding = c(100, 100),
          MahaConstant = 1,
          Method = "Mahal",
          XandY = TRUE,
          Second = TRUE,
          MaxValueRun = 5,
          NAsAroundRun = c(2, 2),
          LogFile = paste0(tempdir(), "/ArtifactCleanupLog.rds")
        ),
        error = function(e2) {
          cat("WARN: clean_artifact Mahal failed, retrying with MAD. Error:", conditionMessage(e2), "\n")
          tryCatch(
            clean_artifact(
              cleanedData,
              MADWindow = 20,
              MADConstant = 4,
              MADPadding = c(100, 100),
              MahaConstant = 1,
              Method = "MAD",
              XandY = TRUE,
              Second = TRUE,
              MaxValueRun = 5,
              NAsAroundRun = c(2, 2),
              LogFile = paste0(tempdir(), "/ArtifactCleanupLog.rds")
            ),
            error = function(e3) {
              cat("WARN: clean_artifact MAD failed, skipping artifact cleanup. Error:", conditionMessage(e3), "\n")
              cleanedData
            }
          )
        }
      )
    }
  )

  # make sure usable samples is large enough
  cleanedData <- rm_sparse_events(
    data = cleanedData,
    BaselineWindow = c(0, 7000),
    CriticalWindow = c(7000, 12000),
    BaselineRequired = 60,
    CriticalRequired = 60
  )

  datlinear <- interpolate_NAs(
    cleanedData,
    Method = "linear",
    XandY = TRUE,
    MinData = 2
  )

  # change units
  reference_pupil_mm <- 8
  reference_pupil_au <- 6971
  scale_factor <- reference_pupil_mm / reference_pupil_au

  datlinear_mm <- datlinear %>%
    mutate(Pupil = Pupil * scale_factor)

  # wavelength
  datfilter <- apply_butter(
    datlinear_mm,
    n = 2,
    W = 0.05,
    type = "low",
    plane = "z"
  )

  # Baseline identification
  datafinal <- baseline(
    datfilter,
    BaselineWindow = c(0, 7000),
    BaselineType = "Subtraction"
  )

  # save data & downsampling
  check_samplingrate(datafinal)
  downsampled_all <- downsample(datafinal, SamplingRate = 1000, NewRate = 25)

  # temporal smoothing
  smooth_window <- 5  
  downsampled_all <- downsampled_all %>%
    arrange(Subject, Event, Time) %>%
    group_by(Subject, Event) %>%
    mutate(
      Pupil_smoothed = as.numeric(stats::filter(
        Pupil,
        rep(1 / smooth_window, smooth_window),
        sides = 2
      )),
      Gaze_X_smoothed = as.numeric(stats::filter(
        Gaze_X,
        rep(1 / smooth_window, smooth_window),
        sides = 2
      )),
      Gaze_Y_smoothed = as.numeric(stats::filter(
        Gaze_Y,
        rep(1 / smooth_window, smooth_window),
        sides = 2
      )),
      Pupil = dplyr::coalesce(Pupil_smoothed, Pupil),
      Gaze_X = dplyr::coalesce(Gaze_X_smoothed, Gaze_X),
      Gaze_Y = dplyr::coalesce(Gaze_Y_smoothed, Gaze_Y)
    ) %>%
    select(-Pupil_smoothed, -Gaze_X_smoothed, -Gaze_Y_smoothed) %>%
    ungroup()

  # save data
  saveRDS(downsampled_all, output_file)
  cat("✅ Saved:", output_file, "\n")
  invisible(output_file)
}

# Wrap single-file preprocessing so the batch run can continue after errors.
safe_process_one <- function(file, output_dir) {
  tryCatch(
    process_one(file, output_dir),
    error = function(e) {
      cat("ERROR:", file, "-", conditionMessage(e), "\n")
      return(NA_character_)
    }
  )
}

results <- lapply(files, safe_process_one, output_dir = output_dir)
cat("Done. Files processed:", length(results), "\n")
