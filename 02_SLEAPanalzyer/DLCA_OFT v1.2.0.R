# ================================================================
# SLEAP/DLC Open Field Test Analyzer
# QC-first OFT analysis with frame-level and summary outputs
# Author: Tobias Pohl / ChatGPT revision
# Version: 1.2.0
# Date: 2026-05-11
# ================================================================
#
# Main changes vs. v1.1.0:
# - Centralized configuration
# - No hard-coded working directory dependency
# - Tracking QC before and after cleaning
# - Short-gap interpolation using DLCAnalyzer CleanTrackingData()
# - Frame-level ethogram output
# - Center/periphery/corner bout metrics
# - Wall-distance/thigmotaxis metrics
# - Time-binned summaries for habituation/exploration trajectories
# - Failure log for files that could not be processed
#
# Biological interpretation note:
# OFT should be read as a profile rather than a single "anxiety" score.
# Always interpret center avoidance together with locomotion, immobility,
# tracking quality, and the time course of exploration.
# ================================================================

# -------------------------------
# 0) Packages
# -------------------------------

required_packages <- c(
  "dplyr", "purrr", "stringr", "readr", "ggplot2", "zoo", "sp", "tibble", "tidyr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# -------------------------------
# 1) User parameters
# -------------------------------

config <- list(
  fps = 30,
  batches = c("B1"),

  # Repository/script setup.
  functions_file = "DLCAnalyzer_Functions_final.R",

  # Experiment paths.
  behavior_root = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior",
  animal_id_code_file = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Planning/animalIDCode.txt",

  # Input/output structure inside each batch folder.
  input_subdir = file.path("OFT", "SLEAP", "formatted"),
  output_subdir = file.path("OFT", "SLEAP", "output_v1.2.0"),
  csv_pattern = "\\.csv$",

  # Filename parsing.
  # Anchored to the first 4 alphanumeric characters to avoid accidental matches later in the filename.
  code_regex = "^[A-Za-z0-9]{4}",

  # Arena and zone geometry.
  arena_size_cm = 49,
  corner_points = c("tl", "tr", "br", "bl"),
  center_scale = 0.5,
  periphery_scale = 0.8,
  corner_scale = 0.4,

  # Tracking cleanup and QC.
  likelihood_cutoff = 0.90,
  max_jump_cm = 15,
  high_missing_fraction = 0.20,
  high_low_likelihood_fraction = 0.20,

  # Movement analysis.
  movement_cutoff_cm_s = 5,
  integration_period_frames = 5,
  immobility_min_duration_s = 1.0,

  # Time-course output.
  bin_size_s = 60,

  # Plots.
  save_plots = TRUE,
  plot_width = 7,
  plot_height = 4
)

# -------------------------------
# 2) Setup
# -------------------------------

get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(path) && nzchar(path)) return(dirname(path))
  }

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  getwd()
}

script_dir <- get_script_dir()
functions_candidates <- c(
  file.path(script_dir, config$functions_file),
  file.path(getwd(), config$functions_file),
  file.path(getwd(), "SLEAPanalyzer", "SLEAPanalzyer", config$functions_file)
)
functions_path <- functions_candidates[file.exists(functions_candidates)][1]

if (is.na(functions_path)) {
  stop(
    "Could not find functions file: ", config$functions_file,
    "\nOpen/source this script from the SLEAPanalzyer directory or update config$functions_file."
  )
}

source(functions_path)

if (file.exists(config$animal_id_code_file)) {
  animalIDCode <- read.table(config$animal_id_code_file, header = TRUE)
} else {
  warning("animal_id_code_file not found: ", config$animal_id_code_file)
  animalIDCode <- data.frame(Code = character(), ID = character())
}

# -------------------------------
# 3) Helper functions
# -------------------------------

recalculate_median_data <- function(tracking) {
  med <- purrr::imap_dfr(
    tracking$data,
    ~ data.frame(
      PointName = .y,
      x = median(.x$x, na.rm = TRUE),
      y = median(.x$y, na.rm = TRUE)
    )
  )
  rownames(med) <- med$PointName
  tracking$median.data <- med
  tracking
}

safe_report <- function(report, name) {
  if (!is.null(report[[name]])) report[[name]] else NA_real_
}

entry_count <- function(x) {
  x <- tidyr::replace_na(as.logical(x), FALSE)
  sum(x & !dplyr::lag(x, default = FALSE), na.rm = TRUE)
}

latency_s <- function(x, fps) {
  x <- tidyr::replace_na(as.logical(x), FALSE)
  idx <- which(x)[1]
  if (is.na(idx)) NA_real_ else (idx - 1) / fps
}

bout_summary <- function(x, fps, min_duration_s = 0) {
  x <- tidyr::replace_na(as.logical(x), FALSE)
  if (length(x) == 0) {
    return(tibble(
      bouts = 0L,
      time_s = 0,
      mean_bout_s = NA_real_,
      max_bout_s = NA_real_,
      latency_s = NA_real_
    ))
  }

  r <- rle(x)
  durations_s <- r$lengths[r$values] / fps
  durations_s <- durations_s[durations_s >= min_duration_s]

  tibble(
    bouts = length(durations_s),
    time_s = sum(durations_s),
    mean_bout_s = ifelse(length(durations_s) == 0, NA_real_, mean(durations_s)),
    max_bout_s = ifelse(length(durations_s) == 0, NA_real_, max(durations_s)),
    latency_s = latency_s(x, fps)
  )
}

qc_tracking <- function(tracking, stage, config) {
  purrr::imap_dfr(tracking$data, function(dat, point) {
    dx <- c(NA_real_, diff(dat$x))
    dy <- c(NA_real_, diff(dat$y))
    step_distance <- sqrt(dx^2 + dy^2)

    tibble(
      stage = stage,
      point = point,
      frames = nrow(dat),
      missing_fraction = mean(is.na(dat$x) | is.na(dat$y)),
      low_likelihood_fraction = if ("likelihood" %in% names(dat)) {
        mean(dat$likelihood < config$likelihood_cutoff, na.rm = TRUE)
      } else {
        NA_real_
      },
      impossible_jump_n = sum(step_distance > config$max_jump_cm, na.rm = TRUE),
      median_likelihood = if ("likelihood" %in% names(dat)) {
        median(dat$likelihood, na.rm = TRUE)
      } else {
        NA_real_
      }
    )
  })
}

validate_tracking <- function(tracking, file_id, config) {
  required <- c(config$corner_points, "bodycentre")
  missing_points <- setdiff(required, names(tracking$data))

  if (length(missing_points) > 0) {
    stop(file_id, " is missing required point(s): ", paste(missing_points, collapse = ", "))
  }

  invisible(TRUE)
}

make_frame_table <- function(tracking, file_id, batch, code, config) {
  body <- tracking$data$bodycentre

  center <- IsInZone(tracking, "bodycentre", "center")
  periphery <- IsInZone(tracking, "bodycentre", "periphery", invert = TRUE)
  corners <- IsInZone(tracking, "bodycentre", tracking$corner.names)
  arena <- IsInZone(tracking, "bodycentre", "arena")
  wall_distance <- GetDistanceToZoneBorder(tracking, zone = "arena", point = "bodycentre")

  tibble(
    file = file_id,
    Batch = batch,
    Code = code,
    frame = body$frame,
    time_s = frame / config$fps,
    x = body$x,
    y = body$y,
    likelihood = body$likelihood,
    speed_cm_s = body$speed * config$fps,
    acceleration_cm_s2 = body$acceleration * config$fps * config$fps,
    is_moving = body$is.moving,
    is_immobile = !body$is.moving,
    in_arena = arena,
    in_center = center,
    in_periphery = periphery,
    in_corners = corners,
    wall_distance_cm = wall_distance
  )
}

make_binned_table <- function(frame_tbl, config) {
  frame_tbl %>%
    mutate(bin_s = floor(time_s / config$bin_size_s) * config$bin_size_s) %>%
    group_by(file, Batch, Code, bin_s) %>%
    summarise(
      bin_start_s = min(time_s, na.rm = TRUE),
      bin_end_s = max(time_s, na.rm = TRUE),
      frames = n(),
      center_time_s = sum(in_center, na.rm = TRUE) / config$fps,
      periphery_time_s = sum(in_periphery, na.rm = TRUE) / config$fps,
      corner_time_s = sum(in_corners, na.rm = TRUE) / config$fps,
      moving_time_s = sum(is_moving, na.rm = TRUE) / config$fps,
      immobile_time_s = sum(is_immobile, na.rm = TRUE) / config$fps,
      distance_cm = sum(speed_cm_s / config$fps, na.rm = TRUE),
      mean_speed_cm_s = mean(speed_cm_s, na.rm = TRUE),
      mean_wall_distance_cm = mean(wall_distance_cm, na.rm = TRUE),
      center_entries = entry_count(in_center),
      .groups = "drop"
    )
}

make_summary_row <- function(tracking, frame_tbl, file_id, batch, code, animalIDCode, qc_tbl, config) {
  total_time_s <- nrow(frame_tbl) / config$fps
  total_distance_cm <- safe_report(tracking$Report, "bodycentre.raw.distance")

  center_bouts <- bout_summary(frame_tbl$in_center, config$fps)
  periphery_bouts <- bout_summary(frame_tbl$in_periphery, config$fps)
  corner_bouts <- bout_summary(frame_tbl$in_corners, config$fps)
  immobile_bouts <- bout_summary(
    frame_tbl$is_immobile,
    config$fps,
    min_duration_s = config$immobility_min_duration_s
  )

  id_match <- animalIDCode[animalIDCode$Code == code, "ID"]
  animal_id <- if (length(id_match) == 0) NA else id_match[1]

  qc_body <- qc_tbl %>% filter(stage == "raw_calibrated", point == "bodycentre")

  tibble(
    file = file_id,
    ID = animal_id,
    Batch = batch,
    Code = code,
    fps = config$fps,
    total_time_s = total_time_s,

    distance_cm = total_distance_cm,
    mean_speed_cm_s = safe_report(tracking$Report, "bodycentre.raw.speed"),
    moving_speed_cm_s = safe_report(tracking$Report, "bodycentre.speed.moving"),
    moving_time_s = safe_report(tracking$Report, "bodycentre.time.moving"),
    stationary_time_s = safe_report(tracking$Report, "bodycentre.time.stationary"),
    percent_moving = safe_report(tracking$Report, "bodycentre.percentage.moving"),

    center_time_s = safe_report(tracking$Report, "bodycentre.center.total.time"),
    center_time_percent = center_time_s / total_time_s * 100,
    center_distance_cm = safe_report(tracking$Report, "bodycentre.center.raw.distance"),
    center_distance_percent = center_distance_cm / total_distance_cm * 100,
    center_entries = entry_count(frame_tbl$in_center),
    center_transitions = safe_report(tracking$Report, "bodycentre.center.transitions"),
    center_latency_s = center_bouts$latency_s,
    center_mean_bout_s = center_bouts$mean_bout_s,
    center_max_bout_s = center_bouts$max_bout_s,

    periphery_time_s = safe_report(tracking$Report, "bodycentre.periphery.total.time"),
    periphery_time_percent = periphery_time_s / total_time_s * 100,
    periphery_distance_cm = safe_report(tracking$Report, "bodycentre.periphery.raw.distance"),
    periphery_entries = entry_count(frame_tbl$in_periphery),
    periphery_mean_bout_s = periphery_bouts$mean_bout_s,
    periphery_max_bout_s = periphery_bouts$max_bout_s,

    corner_time_s = safe_report(tracking$Report, "bodycentre.corners.total.time"),
    corner_time_percent = corner_time_s / total_time_s * 100,
    corner_entries = entry_count(frame_tbl$in_corners),
    corner_mean_bout_s = corner_bouts$mean_bout_s,
    corner_max_bout_s = corner_bouts$max_bout_s,

    immobility_bouts = immobile_bouts$bouts,
    immobility_time_s = immobile_bouts$time_s,
    immobility_mean_bout_s = immobile_bouts$mean_bout_s,
    immobility_max_bout_s = immobile_bouts$max_bout_s,

    mean_wall_distance_cm = mean(frame_tbl$wall_distance_cm, na.rm = TRUE),
    median_wall_distance_cm = median(frame_tbl$wall_distance_cm, na.rm = TRUE),
    min_wall_distance_cm = min(frame_tbl$wall_distance_cm, na.rm = TRUE),
    max_wall_distance_cm = max(frame_tbl$wall_distance_cm, na.rm = TRUE),

    bodycentre_missing_fraction = qc_body$missing_fraction[1],
    bodycentre_low_likelihood_fraction = qc_body$low_likelihood_fraction[1],
    bodycentre_impossible_jump_n = qc_body$impossible_jump_n[1],
    qc_flag = bodycentre_missing_fraction > config$high_missing_fraction |
      bodycentre_low_likelihood_fraction > config$high_low_likelihood_fraction |
      bodycentre_impossible_jump_n > 0
  )
}

save_optional_plots <- function(tracking, file_id, plot_dir, config) {
  if (!isTRUE(config$save_plots)) return(invisible(NULL))

  density_plots <- PlotDensityPaths(tracking, points = c("bodycentre"))
  ggsave(
    filename = file.path(plot_dir, paste0(file_id, "_DensityPath.png")),
    plot = density_plots$bodycentre,
    width = config$plot_width,
    height = config$plot_height
  )

  zone_plot <- PlotZones(tracking)
  ggsave(
    filename = file.path(plot_dir, paste0(file_id, "_Zones.png")),
    plot = zone_plot,
    width = config$plot_width,
    height = config$plot_height
  )

  invisible(NULL)
}

process_oft_file <- function(file, batch, dirs, animalIDCode, config) {
  input_file <- file.path(dirs$input, file)
  file_id <- sub("\\.csv$", "", basename(input_file))
  code <- stringr::str_extract(file_id, config$code_regex)

  tracking <- ReadDLCDataFromCSV(file = input_file, fps = config$fps)
  validate_tracking(tracking, file_id, config)

  tracking <- CalibrateTrackingData(
    tracking,
    method = "area",
    in.metric = config$arena_size_cm * config$arena_size_cm,
    points = config$corner_points
  )

  tracking <- recalculate_median_data(tracking)

  tracking <- AddOFTZones(
    tracking,
    points = config$corner_points,
    scale_center = config$center_scale,
    scale_periphery = config$periphery_scale,
    scale_corners = config$corner_scale
  )

  qc_raw <- qc_tracking(tracking, "raw_calibrated", config)

  tracking <- CleanTrackingData(
    tracking,
    likelihoodcutoff = config$likelihood_cutoff,
    existence.pol = tracking$zones$arena,
    maxdelta = config$max_jump_cm
  )

  tracking <- recalculate_median_data(tracking)

  tracking <- OFTAnalysis(
    tracking,
    points = "bodycentre",
    movement_cutoff = config$movement_cutoff_cm_s,
    integration_period = config$integration_period_frames
  )

  qc_clean <- qc_tracking(tracking, "cleaned", config)
  qc_tbl <- bind_rows(qc_raw, qc_clean)
  frame_tbl <- make_frame_table(tracking, file_id, batch, code, config)
  binned_tbl <- make_binned_table(frame_tbl, config)
  summary_tbl <- make_summary_row(tracking, frame_tbl, file_id, batch, code, animalIDCode, qc_tbl, config)

  readr::write_csv(summary_tbl, file.path(dirs$tables, paste0(file_id, "_summary.csv")))
  readr::write_csv(frame_tbl, file.path(dirs$frames, paste0(file_id, "_frames.csv")))
  readr::write_csv(binned_tbl, file.path(dirs$bins, paste0(file_id, "_bins.csv")))
  readr::write_csv(qc_tbl, file.path(dirs$qc, paste0(file_id, "_qc.csv")))

  save_optional_plots(tracking, file_id, dirs$plots, config)

  list(
    summary = summary_tbl,
    frames = frame_tbl,
    bins = binned_tbl,
    qc = qc_tbl
  )
}

# -------------------------------
# 4) Batch processing
# -------------------------------

all_summaries <- list()
all_bins <- list()
all_qc <- list()
failures <- list()

for (batch in config$batches) {
  dirs <- list(
    input = file.path(config$behavior_root, batch, config$input_subdir),
    output = file.path(config$behavior_root, batch, config$output_subdir)
  )

  dirs$tables <- file.path(dirs$output, "tables")
  dirs$frames <- file.path(dirs$output, "frames")
  dirs$bins <- file.path(dirs$output, "bins")
  dirs$qc <- file.path(dirs$output, "qc")
  dirs$plots <- file.path(dirs$output, "plots")

  purrr::walk(dirs[-1], ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

  file_list <- list.files(path = dirs$input, pattern = config$csv_pattern)

  if (length(file_list) == 0) {
    warning("No CSV files found for batch ", batch, " in ", dirs$input)
    next
  }

  for (file in file_list) {
    message("Processing ", batch, " / ", file)

    result <- tryCatch(
      process_oft_file(file, batch, dirs, animalIDCode, config),
      error = function(e) {
        warning("Failed to process ", file, ": ", conditionMessage(e))
        failures[[length(failures) + 1]] <<- tibble(
          Batch = batch,
          file = file,
          error = conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(result)) {
      all_summaries[[length(all_summaries) + 1]] <- result$summary
      all_bins[[length(all_bins) + 1]] <- result$bins
      all_qc[[length(all_qc) + 1]] <- result$qc %>%
        mutate(Batch = batch, file = sub("\\.csv$", "", file), .before = 1)
    }
  }

  batch_summaries <- bind_rows(all_summaries) %>% filter(Batch == batch)
  batch_bins <- bind_rows(all_bins) %>% filter(Batch == batch)
  batch_qc <- bind_rows(all_qc) %>% filter(Batch == batch)

  if (nrow(batch_summaries) > 0) {
    readr::write_csv(batch_summaries, file.path(dirs$output, "combined_summary.csv"))
    readr::write_csv(batch_bins, file.path(dirs$output, "combined_bins.csv"))
    readr::write_csv(batch_qc, file.path(dirs$output, "combined_qc.csv"))
  }

  message("Processing complete for batch ", batch)
}

if (length(all_summaries) > 0) {
  global_output <- file.path(config$behavior_root, "OFT_SLEAP_output_v1.2.0")
  dir.create(global_output, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(bind_rows(all_summaries), file.path(global_output, "all_batches_summary.csv"))
  readr::write_csv(bind_rows(all_bins), file.path(global_output, "all_batches_bins.csv"))
  readr::write_csv(bind_rows(all_qc), file.path(global_output, "all_batches_qc.csv"))
}

if (length(failures) > 0) {
  failure_tbl <- bind_rows(failures)
  failure_output <- file.path(config$behavior_root, "OFT_SLEAP_output_v1.2.0")
  dir.create(failure_output, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(failure_tbl, file.path(failure_output, "failures.csv"))
  warning(nrow(failure_tbl), " file(s) failed. See failures.csv.")
}

message("OFT analysis complete.")
