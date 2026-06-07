# ================================================================
# SLEAP/DLC Social Interaction Analyzer
# Event-based animal-animal interaction analysis with QC
# Author: Tobias Pohl / ChatGPT revision
# ================================================================
#
# What this script adds vs. the exploratory v0.0.2 script:
# - No hard-coded per-file overwrite inside the batch loop
# - Centralized parameters for fps and thresholds
# - Robust interpolation with QC reporting
# - Geometry helpers for distances, headings, and angles
# - Frame-level ethogram / event table
# - Directional interaction metrics: animal1 -> animal2 and animal2 -> animal1
# - Bout statistics: duration, count, latency, mean/max bout duration
# - Proximity, facing, parallel, anti-parallel, following, and avoidance metrics
# - Combined per-file summary table
# - QC plots and distance/event plots
#
# Expected tracking object:
# This script assumes that DLCAnalyzer_Functions_final.R provides:
#   ReadDLCDataFromCSV(file, fps)
# and returns Tracking$data[[bodypart]]$x and Tracking$data[[bodypart]]$y
#
# Bodypart names expected by default:
# animal 1: nose_1, leftEar_1, rightEar_1, bodycentre_1, leftSide_1, rightSide_1, tailBase_1, tailEnd_1
# animal 2: nose_2, leftEar_2, rightEar_2, bodycentre_2, leftSide_2, rightSide_2, tailBase_2, tailEnd_2
#
# ================================================================

# -------------------------------
# 0) Packages
# -------------------------------

required_packages <- c(
  "dplyr", "tidyr", "purrr", "stringr", "readr", "ggplot2",
  "zoo", "tibble", "fs", "rlang"
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
  overwrite_outputs = TRUE,
  recursive_input_search = FALSE,
  min_track_duration_s = 60,
  strict_qc_stop = FALSE,
  min_event_duration_frames = 3,

  # Set these paths.
  # Use forward slashes also on Windows to avoid escaping problems.
  working_dir = "C:/Users/topohl/Documents/GitHub/SLEAPanalyzer/02_SLEAPanalzyer",
  functions_file = "DLCAnalyzer_Functions_final.R",

  input_dir  = "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC",
  output_dir = "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC/output_event_based",

  # File pattern.
  csv_pattern = "\\.csv$",

  # Coordinate units and optional arena calibration.
  # If use_arena_calibration = TRUE, x/y coordinates are converted from pixels to cm
  # before distances, speeds, jumps, and event thresholds are computed.
  unit = "px",
  use_arena_calibration = TRUE,
  arena_geom_dir = "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC/geom",
  arena_geom_suffix = "_locs.csv",
  arena_width_cm = 49,
  arena_height_cm = 49,
  arena_corner_names = c("tl", "tr", "br", "bl"),

  # Interaction thresholds.
  # If use_arena_calibration = TRUE, these thresholds are interpreted as cm.
  # Otherwise, they are interpreted in raw coordinate units, usually pixels.
  contact_dist = 30,
  side_by_side_dist = 80,
  close_proximity_dist = 80,
  medium_proximity_dist = 160,
  proximity_bins = c(0, 5, 10, 20, 40, Inf),
  proximity_bin_labels = c("0_5", "5_10", "10_20", "20_40", "gt_40"),

  # Orientation thresholds.
  # angle_to_target is 0 when the animal nose/body axis points directly to target.
  # Use <= facing_angle for directed investigation.
  facing_angle_deg = 60,
  parallel_angle_deg = 45,
  anti_parallel_angle_deg = 135,

  # Movement thresholds.
  movement_cutoff = 5,
  follow_lag_s = 1.0,
  follow_dist = 120,
  avoidance_lag_s = 1.0,
  avoidance_start_dist = 80,
  avoidance_delta_dist = 40,
  approach_lag_s = 1.0,
  approach_start_dist = 20,
  approach_delta_dist = 5,
  retreat_lag_s = 1.0,
  retreat_start_dist = 10,
  retreat_delta_dist = 5,

  # QC thresholds.
  max_interpolation_gap_frames = 30,
  impossible_jump_dist = 100,
  high_missing_fraction = 0.20,

  # Plots.
  save_plots = TRUE,
  plot_width = 8,
  plot_height = 5
)

# -------------------------------
# 2) Setup
# -------------------------------

script_dir <- if (dir.exists(config$working_dir)) normalizePath(config$working_dir) else getwd()

functions_candidates <- c(
  config$functions_file,
  file.path(script_dir, config$functions_file),
  file.path(getwd(), config$functions_file)
)
functions_path <- functions_candidates[file.exists(functions_candidates)][1]

if (is.na(functions_path)) {
  stop(
    "Could not find functions file: ", config$functions_file,
    "\nChecked: ", paste(functions_candidates, collapse = "; "),
    "\nExpected function ReadDLCDataFromCSV(file, fps)."
  )
}

source(functions_path)

fs::dir_create(config$output_dir)
fs::dir_create(file.path(config$output_dir, "tables"))
fs::dir_create(file.path(config$output_dir, "frames"))
fs::dir_create(file.path(config$output_dir, "qc"))
fs::dir_create(file.path(config$output_dir, "figures"))

# -------------------------------
# 3) Bodypart definitions
# -------------------------------

bodyparts_animal1 <- c(
  "nose_1", "leftEar_1", "rightEar_1", "bodycentre_1",
  "leftSide_1", "rightSide_1", "tailBase_1", "tailEnd_1"
)

bodyparts_animal2 <- c(
  "nose_2", "leftEar_2", "rightEar_2", "bodycentre_2",
  "leftSide_2", "rightSide_2", "tailBase_2", "tailEnd_2"
)

bodyparts_all <- c(bodyparts_animal1, bodyparts_animal2)

# -------------------------------
# 4) Utility functions
# -------------------------------

safe_divide <- function(x, y) {
  out <- rep(NA_real_, length.out = max(length(x), length(y)))
  ok <- !is.na(y) & y != 0
  out[ok] <- x[ok] / y[ok]
  out
}

suppress_short_events <- function(event_vec, min_frames = 1) {
  event_vec <- tidyr::replace_na(as.logical(event_vec), FALSE)
  if (min_frames <= 1 || length(event_vec) == 0) return(event_vec)

  r <- rle(event_vec)
  keep <- r$values & r$lengths >= min_frames
  inverse.rle(list(lengths = r$lengths, values = keep))
}

vec_norm <- function(x, y) {
  sqrt(x^2 + y^2)
}

angle_between_vectors_deg <- function(ax, ay, bx, by) {
  dot <- ax * bx + ay * by
  mag_a <- vec_norm(ax, ay)
  mag_b <- vec_norm(bx, by)

  cos_theta <- safe_divide(dot, mag_a * mag_b)
  cos_theta <- pmin(1, pmax(-1, cos_theta))

  angle <- acos(cos_theta) * 180 / pi
  angle[is.nan(angle)] <- NA_real_
  angle
}

circular_angle_diff_deg <- function(a, b) {
  # Absolute smallest difference between two angles in degrees.
  diff <- abs(a - b) %% 360
  pmin(diff, 360 - diff)
}

get_xy <- function(Tracking, bodypart) {
  if (!bodypart %in% names(Tracking$data)) {
    stop("Bodypart not found in Tracking$data: ", bodypart)
  }

  tibble(
    x = as.numeric(Tracking$data[[bodypart]]$x),
    y = as.numeric(Tracking$data[[bodypart]]$y)
  )
}

# ---- Arena calibration helpers ----
set_xy <- function(Tracking, bodypart, x, y) {
  Tracking$data[[bodypart]]$x <- x
  Tracking$data[[bodypart]]$y <- y
  Tracking
}

get_output_unit <- function(config) {
  ifelse(isTRUE(config$use_arena_calibration), "cm", config$unit)
}

find_geom_file <- function(file_id, config) {
  if (is.null(config$arena_geom_dir) || !dir.exists(config$arena_geom_dir)) {
    stop("arena_geom_dir does not exist: ", config$arena_geom_dir)
  }

  candidates <- c(
    file.path(config$arena_geom_dir, paste0(file_id, config$arena_geom_suffix)),
    file.path(config$arena_geom_dir, paste0(file_id, "_geom", config$arena_geom_suffix)),
    file.path(config$arena_geom_dir, paste0(file_id, "_arena", config$arena_geom_suffix))
  )

  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) > 0) return(candidates[1])

  geom_files <- list.files(config$arena_geom_dir, pattern = "\\.csv$", full.names = TRUE)
  loose <- geom_files[stringr::str_detect(basename(geom_files), stringr::fixed(file_id))]
  if (length(loose) > 0) return(loose[1])

  stop("No arena geometry CSV found for ", file_id, " in ", config$arena_geom_dir)
}

read_static_arena_geometry <- function(file_id, config) {
  geom_file <- find_geom_file(file_id, config)
  geom_tbl <- readr::read_csv(geom_file, show_col_types = FALSE)

  required_cols <- as.vector(t(outer(config$arena_corner_names, c("x", "y"), paste, sep = "_")))
  missing_cols <- setdiff(required_cols, names(geom_tbl))
  if (length(missing_cols) > 0) {
    stop("Arena geometry file is missing columns: ", paste(missing_cols, collapse = ", "), " | file: ", geom_file)
  }

  purrr::map_dfr(config$arena_corner_names, function(corner) {
    tibble(
      corner = corner,
      x = median(geom_tbl[[paste0(corner, "_x")]], na.rm = TRUE),
      y = median(geom_tbl[[paste0(corner, "_y")]], na.rm = TRUE),
      geom_file = geom_file
    )
  })
}

estimate_arena_calibration <- function(arena_geom, config) {
  get_corner <- function(name) arena_geom |> filter(corner == name) |> slice(1)

  tl <- get_corner("tl")
  tr <- get_corner("tr")
  br <- get_corner("br")
  bl <- get_corner("bl")

  top_px <- vec_norm(tr$x - tl$x, tr$y - tl$y)
  bottom_px <- vec_norm(br$x - bl$x, br$y - bl$y)
  left_px <- vec_norm(bl$x - tl$x, bl$y - tl$y)
  right_px <- vec_norm(br$x - tr$x, br$y - tr$y)

  width_px <- mean(c(top_px, bottom_px), na.rm = TRUE)
  height_px <- mean(c(left_px, right_px), na.rm = TRUE)

  cm_per_px_x <- safe_divide(config$arena_width_cm, width_px)
  cm_per_px_y <- safe_divide(config$arena_height_cm, height_px)
  cm_per_px_mean <- mean(c(cm_per_px_x, cm_per_px_y), na.rm = TRUE)

  if (is.na(cm_per_px_mean) || cm_per_px_mean <= 0) {
    stop("Invalid arena calibration. Check arena corner coordinates.")
  }

  tibble(
    width_px = width_px,
    height_px = height_px,
    arena_width_cm = config$arena_width_cm,
    arena_height_cm = config$arena_height_cm,
    cm_per_px_x = cm_per_px_x,
    cm_per_px_y = cm_per_px_y,
    cm_per_px_mean = cm_per_px_mean,
    calibration_anisotropy_percent = abs(cm_per_px_x - cm_per_px_y) / cm_per_px_mean * 100
  )
}

calibrate_tracking_to_cm <- function(Tracking, bodyparts, arena_geom, calibration) {
  origin_x <- median(arena_geom$x, na.rm = TRUE)
  origin_y <- median(arena_geom$y, na.rm = TRUE)

  for (bp in bodyparts) {
    if (!bp %in% names(Tracking$data)) next

    x_px <- as.numeric(Tracking$data[[bp]]$x)
    y_px <- as.numeric(Tracking$data[[bp]]$y)

    Tracking <- set_xy(
      Tracking,
      bp,
      x = (x_px - origin_x) * calibration$cm_per_px_x[1],
      y = (y_px - origin_y) * calibration$cm_per_px_y[1]
    )
  }

  Tracking
}

maybe_calibrate_tracking <- function(Tracking, file_id, config, bodyparts_all) {
  if (!isTRUE(config$use_arena_calibration)) {
    return(list(
      Tracking = Tracking,
      arena_geom = tibble(),
      calibration = tibble(
        file = file_id,
        calibrated = FALSE,
        unit_after_calibration = config$unit
      )
    ))
  }

  arena_geom <- read_static_arena_geometry(file_id, config)
  calibration <- estimate_arena_calibration(arena_geom, config) |>
    mutate(
      file = file_id,
      calibrated = TRUE,
      unit_after_calibration = "cm",
      geom_file = unique(arena_geom$geom_file)[1],
      .before = 1
    )

  Tracking <- calibrate_tracking_to_cm(Tracking, bodyparts_all, arena_geom, calibration)

  list(
    Tracking = Tracking,
    arena_geom = arena_geom |> mutate(file = file_id, .before = 1),
    calibration = calibration
  )
}

distance_bp <- function(Tracking, bp1, bp2) {
  a <- get_xy(Tracking, bp1)
  b <- get_xy(Tracking, bp2)
  vec_norm(a$x - b$x, a$y - b$y)
}

heading_vector <- function(Tracking, animal_id) {
  nose <- get_xy(Tracking, paste0("nose_", animal_id))
  body <- get_xy(Tracking, paste0("bodycentre_", animal_id))

  tibble(
    x = nose$x - body$x,
    y = nose$y - body$y
  )
}

target_vector <- function(Tracking, from_bp, to_bp) {
  from <- get_xy(Tracking, from_bp)
  to   <- get_xy(Tracking, to_bp)

  tibble(
    x = to$x - from$x,
    y = to$y - from$y
  )
}

angle_animal_to_bodypart <- function(Tracking, animal_id, target_bp) {
  hv <- heading_vector(Tracking, animal_id)
  tv <- target_vector(Tracking, paste0("nose_", animal_id), target_bp)

  angle_between_vectors_deg(hv$x, hv$y, tv$x, tv$y)
}

heading_angle_deg <- function(Tracking, animal_id) {
  hv <- heading_vector(Tracking, animal_id)
  atan2(hv$y, hv$x) * 180 / pi
}

speed_bp <- function(Tracking, bp, fps) {
  xy <- get_xy(Tracking, bp)
  dx <- c(NA_real_, diff(xy$x))
  dy <- c(NA_real_, diff(xy$y))
  vec_norm(dx, dy) * fps
}

first_true_time <- function(event_vec, fps) {
  idx <- which(event_vec %in% TRUE)[1]
  if (is.na(idx)) return(NA_real_)
  (idx - 1) / fps
}

summarise_event <- function(event_vec, fps, min_event_duration_frames = 1) {
  event_vec <- suppress_short_events(event_vec, min_frames = min_event_duration_frames)

  starts <- event_vec & !dplyr::lag(event_vec, default = FALSE)

  bout_id <- cumsum(starts)
  bout_id[!event_vec] <- NA_integer_

  bout_df <- tibble(active = event_vec, bout_id = bout_id) |>
    filter(active) |>
    group_by(bout_id) |>
    summarise(frames = n(), .groups = "drop")

  tibble(
    duration_s = sum(event_vec, na.rm = TRUE) / fps,
    percent_time = 100 * mean(event_vec, na.rm = TRUE),
    bout_n = sum(starts, na.rm = TRUE),
    latency_s = first_true_time(event_vec, fps),
    mean_bout_s = ifelse(nrow(bout_df) == 0, NA_real_, mean(bout_df$frames) / fps),
    max_bout_s  = ifelse(nrow(bout_df) == 0, NA_real_, max(bout_df$frames) / fps),
    event_frames = sum(event_vec, na.rm = TRUE)
  )
}

interbout_intervals_s <- function(event_vec, fps, min_event_duration_frames = 1) {
  event_vec <- suppress_short_events(event_vec, min_frames = min_event_duration_frames)
  if (!any(event_vec, na.rm = TRUE)) return(numeric(0))

  starts <- which(event_vec & !dplyr::lag(event_vec, default = FALSE))
  ends <- which(event_vec & !dplyr::lead(event_vec, default = FALSE))

  if (length(starts) <= 1 || length(ends) == 0) return(numeric(0))

  intervals_frames <- starts[-1] - ends[-length(ends)] - 1
  intervals_frames <- intervals_frames[intervals_frames >= 0]
  intervals_frames / fps
}

interpolate_with_qc <- function(x, maxgap = Inf) {
  x_orig <- as.numeric(x)
  n <- length(x_orig)

  # All values missing: return unchanged and let QC flag it.
  if (all(is.na(x_orig))) {
    return(list(
      x = x_orig,
      n_missing_original = n,
      frac_missing_original = 1,
      n_interpolated = 0,
      frac_interpolated = 0,
      longest_missing_run = n
    ))
  }

  missing_orig <- is.na(x_orig)
  r <- rle(missing_orig)
  longest_missing_run <- ifelse(any(r$values), max(r$lengths[r$values]), 0)

  # Fill leading/trailing NA with nearest observed value.
  first_non_na <- which(!is.na(x_orig))[1]
  last_non_na  <- tail(which(!is.na(x_orig)), 1)

  x_filled <- x_orig
  if (first_non_na > 1) x_filled[1:(first_non_na - 1)] <- x_orig[first_non_na]
  if (last_non_na < n) x_filled[(last_non_na + 1):n] <- x_orig[last_non_na]

  # Interpolate internal gaps up to maxgap.
  x_interp <- zoo::na.approx(x_filled, na.rm = FALSE, maxgap = maxgap)

  interpolated <- missing_orig & !is.na(x_interp)

  list(
    x = x_interp,
    n_missing_original = sum(missing_orig),
    frac_missing_original = mean(missing_orig),
    n_interpolated = sum(interpolated),
    frac_interpolated = mean(interpolated),
    longest_missing_run = longest_missing_run
  )
}

qc_tracking <- function(Tracking, bodyparts, config, file_id) {
  qc_rows <- list()

  for (bp in bodyparts) {
    if (!bp %in% names(Tracking$data)) {
      qc_rows[[length(qc_rows) + 1]] <- tibble(
        file = file_id,
        bodypart = bp,
        coordinate = NA_character_,
        n_missing_original = NA_integer_,
        frac_missing_original = NA_real_,
        n_interpolated = NA_integer_,
        frac_interpolated = NA_real_,
        longest_missing_run = NA_integer_,
        impossible_jump_n = NA_integer_,
        flag_high_missing = TRUE,
        flag_missing_bodypart = TRUE,
        qc_fail = TRUE
      )
      next
    }

    for (coord in c("x", "y")) {
      raw <- as.numeric(Tracking$data[[bp]][[coord]])
      jump_n <- sum(abs(c(NA_real_, diff(raw))) > config$impossible_jump_dist, na.rm = TRUE)

      interp <- interpolate_with_qc(raw, maxgap = config$max_interpolation_gap_frames)
      Tracking$data[[bp]][[coord]] <- interp$x

      qc_rows[[length(qc_rows) + 1]] <- tibble(
        file = file_id,
        bodypart = bp,
        coordinate = coord,
        n_missing_original = interp$n_missing_original,
        frac_missing_original = interp$frac_missing_original,
        n_interpolated = interp$n_interpolated,
        frac_interpolated = interp$frac_interpolated,
        longest_missing_run = interp$longest_missing_run,
        impossible_jump_n = jump_n,
        flag_high_missing = interp$frac_missing_original > config$high_missing_fraction,
        flag_missing_bodypart = FALSE,
        qc_fail = interp$frac_missing_original > config$high_missing_fraction || jump_n > 0
      )
    }
  }

  list(
    Tracking = Tracking,
    qc = bind_rows(qc_rows)
  )
}

validate_tracking <- function(Tracking, required_bodyparts, config, file_id) {
  missing_bodyparts <- setdiff(required_bodyparts, names(Tracking$data))
  if (length(missing_bodyparts) > 0) {
    stop("Missing required bodyparts in ", file_id, ": ", paste(missing_bodyparts, collapse = ", "))
  }

  n_frames <- length(Tracking$data[[required_bodyparts[1]]]$x)
  duration_s <- n_frames / config$fps
  if (is.na(duration_s) || duration_s < config$min_track_duration_s) {
    stop(file_id, " is too short for reliable social-interaction analysis: ", round(duration_s, 2), " s")
  }

  invisible(TRUE)
}

validate_qc <- function(qc_tbl, config, file_id) {
  essential <- qc_tbl |>
    filter(bodypart %in% c("nose_1", "nose_2", "bodycentre_1", "bodycentre_2"))

  hard_fail <- any(essential$qc_fail, na.rm = TRUE)
  if (isTRUE(config$strict_qc_stop) && hard_fail) {
    stop(file_id, " failed strict QC for one or more essential bodyparts.")
  }

  invisible(TRUE)
}

make_event_table <- function(Tracking, config, file_id) {
  fps <- config$fps
  n_frames <- length(Tracking$data$nose_1$x)
  required_lengths <- purrr::map_int(c("nose_1", "nose_2", "bodycentre_1", "bodycentre_2"), ~ length(Tracking$data[[.x]]$x))
  if (length(unique(required_lengths)) != 1) {
    stop(file_id, " has inconsistent frame counts across required bodyparts.")
  }

  frame_tbl <- tibble(
    file = file_id,
    frame = seq_len(n_frames),
    time_s = (frame - 1) / fps
  )

  # Core distances.
  d_nose_nose <- distance_bp(Tracking, "nose_1", "nose_2")
  d_body_body <- distance_bp(Tracking, "bodycentre_1", "bodycentre_2")
  d_tail_tail <- distance_bp(Tracking, "tailBase_1", "tailBase_2")

  d_nose1_body2 <- distance_bp(Tracking, "nose_1", "bodycentre_2")
  d_nose1_tail2 <- distance_bp(Tracking, "nose_1", "tailBase_2")
  d_nose2_body1 <- distance_bp(Tracking, "nose_2", "bodycentre_1")
  d_nose2_tail1 <- distance_bp(Tracking, "nose_2", "tailBase_1")

  d_tail1_nose2 <- distance_bp(Tracking, "tailBase_1", "nose_2")
  d_tail2_nose1 <- distance_bp(Tracking, "tailBase_2", "nose_1")

  # Directional angles.
  a1_to_nose2 <- angle_animal_to_bodypart(Tracking, 1, "nose_2")
  a1_to_body2 <- angle_animal_to_bodypart(Tracking, 1, "bodycentre_2")
  a1_to_tail2 <- angle_animal_to_bodypart(Tracking, 1, "tailBase_2")

  a2_to_nose1 <- angle_animal_to_bodypart(Tracking, 2, "nose_1")
  a2_to_body1 <- angle_animal_to_bodypart(Tracking, 2, "bodycentre_1")
  a2_to_tail1 <- angle_animal_to_bodypart(Tracking, 2, "tailBase_1")

  heading1 <- heading_angle_deg(Tracking, 1)
  heading2 <- heading_angle_deg(Tracking, 2)
  heading_diff <- circular_angle_diff_deg(heading1, heading2)

  speed1 <- speed_bp(Tracking, "bodycentre_1", fps)
  speed2 <- speed_bp(Tracking, "bodycentre_2", fps)

  follow_lag_frames <- max(1, round(config$follow_lag_s * fps))
  avoidance_lag_frames <- max(1, round(config$avoidance_lag_s * fps))
  approach_lag_frames <- max(1, round(config$approach_lag_s * fps))
  retreat_lag_frames <- max(1, round(config$retreat_lag_s * fps))

  # Basic directional contacts.
  a1_nose_nose <- d_nose_nose <= config$contact_dist &
    a1_to_nose2 <= config$facing_angle_deg

  a1_nose_body <- d_nose1_body2 <= config$contact_dist &
    a1_to_body2 <= config$facing_angle_deg

  a1_nose_tail <- d_nose1_tail2 <= config$contact_dist &
    a1_to_tail2 <= config$facing_angle_deg

  a2_nose_nose <- d_nose_nose <= config$contact_dist &
    a2_to_nose1 <= config$facing_angle_deg

  a2_nose_body <- d_nose2_body1 <= config$contact_dist &
    a2_to_body1 <= config$facing_angle_deg

  a2_nose_tail <- d_nose2_tail1 <= config$contact_dist &
    a2_to_tail1 <= config$facing_angle_deg

  # Non-directional posture / proximity.
  side_by_side_parallel <- d_nose_nose <= config$side_by_side_dist &
    d_body_body <= config$side_by_side_dist &
    d_tail_tail <= config$side_by_side_dist &
    heading_diff <= config$parallel_angle_deg

  side_by_side_antiparallel <- d_nose1_tail2 <= config$side_by_side_dist &
    d_body_body <= config$side_by_side_dist &
    d_tail1_nose2 <= config$side_by_side_dist &
    heading_diff >= config$anti_parallel_angle_deg

  close_proximity <- d_body_body <= config$close_proximity_dist
  medium_proximity <- d_body_body > config$close_proximity_dist &
    d_body_body <= config$medium_proximity_dist

  mutually_facing <- a1_to_nose2 <= config$facing_angle_deg &
    a2_to_nose1 <= config$facing_angle_deg &
    d_body_body <= config$medium_proximity_dist

  # Following / avoidance.
  # Interpretation: A1 follows if A2 was moving shortly before, A1 is moving now,
  # and they remain spatially close. This is heuristic and should be validated visually.
  a1_following_a2 <- d_body_body <= config$follow_dist &
    dplyr::lag(speed2, n = follow_lag_frames, default = NA_real_) > config$movement_cutoff &
    speed1 > config$movement_cutoff

  a2_following_a1 <- d_body_body <= config$follow_dist &
    dplyr::lag(speed1, n = follow_lag_frames, default = NA_real_) > config$movement_cutoff &
    speed2 > config$movement_cutoff

  # Avoidance: close now and distance increases after lag.
  future_d <- dplyr::lead(d_body_body, n = avoidance_lag_frames, default = NA_real_)

  a1_avoidance_from_a2 <- d_body_body <= config$avoidance_start_dist &
    future_d - d_body_body >= config$avoidance_delta_dist &
    speed1 > config$movement_cutoff

  a2_avoidance_from_a1 <- d_body_body <= config$avoidance_start_dist &
    future_d - d_body_body >= config$avoidance_delta_dist &
    speed2 > config$movement_cutoff

  # Approach / retreat events based on change in body-centre distance over a short lag.
  past_d_approach <- dplyr::lag(d_body_body, n = approach_lag_frames, default = NA_real_)
  future_d_retreat <- dplyr::lead(d_body_body, n = retreat_lag_frames, default = NA_real_)

  approach_event <- past_d_approach >= config$approach_start_dist &
    past_d_approach - d_body_body >= config$approach_delta_dist &
    (speed1 > config$movement_cutoff | speed2 > config$movement_cutoff)

  retreat_event <- d_body_body <= config$retreat_start_dist &
    future_d_retreat - d_body_body >= config$retreat_delta_dist &
    (speed1 > config$movement_cutoff | speed2 > config$movement_cutoff)

  # Export wide frame-level table.
  frame_tbl <- frame_tbl |>
    mutate(
      d_nose_nose = d_nose_nose,
      d_body_body = d_body_body,
      d_nose1_body2 = d_nose1_body2,
      d_nose1_tail2 = d_nose1_tail2,
      d_nose2_body1 = d_nose2_body1,
      d_nose2_tail1 = d_nose2_tail1,
      angle_a1_to_nose2 = a1_to_nose2,
      angle_a1_to_body2 = a1_to_body2,
      angle_a1_to_tail2 = a1_to_tail2,
      angle_a2_to_nose1 = a2_to_nose1,
      angle_a2_to_body1 = a2_to_body1,
      angle_a2_to_tail1 = a2_to_tail1,
      heading_diff = heading_diff,
      speed1 = speed1,
      speed2 = speed2,

      a1_nose_to_nose2 = suppress_short_events(a1_nose_nose, config$min_event_duration_frames),
      a1_nose_to_body2 = suppress_short_events(a1_nose_body, config$min_event_duration_frames),
      a1_nose_to_tail2 = suppress_short_events(a1_nose_tail, config$min_event_duration_frames),
      a2_nose_to_nose1 = suppress_short_events(a2_nose_nose, config$min_event_duration_frames),
      a2_nose_to_body1 = suppress_short_events(a2_nose_body, config$min_event_duration_frames),
      a2_nose_to_tail1 = suppress_short_events(a2_nose_tail, config$min_event_duration_frames),

      side_by_side_parallel = suppress_short_events(side_by_side_parallel, config$min_event_duration_frames),
      side_by_side_antiparallel = suppress_short_events(side_by_side_antiparallel, config$min_event_duration_frames),
      close_proximity = suppress_short_events(close_proximity, config$min_event_duration_frames),
      medium_proximity = suppress_short_events(medium_proximity, config$min_event_duration_frames),
      mutually_facing = suppress_short_events(mutually_facing, config$min_event_duration_frames),
      a1_following_a2 = suppress_short_events(a1_following_a2, config$min_event_duration_frames),
      a2_following_a1 = suppress_short_events(a2_following_a1, config$min_event_duration_frames),
      a1_avoidance_from_a2 = suppress_short_events(a1_avoidance_from_a2, config$min_event_duration_frames),
      a2_avoidance_from_a1 = suppress_short_events(a2_avoidance_from_a1, config$min_event_duration_frames),
      approach_event = suppress_short_events(approach_event, config$min_event_duration_frames),
      retreat_event = suppress_short_events(retreat_event, config$min_event_duration_frames),

      # Combined directional investigation states.
      a1_face_investigation = a1_nose_to_nose2,
      a2_face_investigation = a2_nose_to_nose1,
      face_investigation = a1_face_investigation | a2_face_investigation,
      a1_body_investigation = a1_nose_to_body2,
      a2_body_investigation = a2_nose_to_body1,
      body_investigation = a1_body_investigation | a2_body_investigation,
      a1_anogenital_investigation = a1_nose_to_tail2,
      a2_anogenital_investigation = a2_nose_to_tail1,
      anogenital_investigation = a1_anogenital_investigation | a2_anogenital_investigation,
      a1_any_directed_contact = a1_face_investigation | a1_body_investigation | a1_anogenital_investigation,
      a2_any_directed_contact = a2_face_investigation | a2_body_investigation | a2_anogenital_investigation,
      mutual_directed_contact = a1_any_directed_contact & a2_any_directed_contact,
      a1_only_directed_contact = a1_any_directed_contact & !a2_any_directed_contact,
      a2_only_directed_contact = a2_any_directed_contact & !a1_any_directed_contact,
      any_social_contact = a1_any_directed_contact | a2_any_directed_contact |
        side_by_side_parallel | side_by_side_antiparallel
    )

  frame_tbl
}

summarise_events_long <- function(frame_tbl, config) {
  fps <- config$fps

  event_cols <- c(
    "a1_nose_to_nose2",
    "a1_nose_to_body2",
    "a1_nose_to_tail2",
    "a2_nose_to_nose1",
    "a2_nose_to_body1",
    "a2_nose_to_tail1",
    "a1_face_investigation",
    "a2_face_investigation",
    "face_investigation",
    "a1_body_investigation",
    "a2_body_investigation",
    "body_investigation",
    "a1_anogenital_investigation",
    "a2_anogenital_investigation",
    "anogenital_investigation",
    "side_by_side_parallel",
    "side_by_side_antiparallel",
    "close_proximity",
    "medium_proximity",
    "mutually_facing",
    "a1_following_a2",
    "a2_following_a1",
    "a1_avoidance_from_a2",
    "a2_avoidance_from_a1",
    "approach_event",
    "retreat_event",
    "a1_any_directed_contact",
    "a2_any_directed_contact",
    "mutual_directed_contact",
    "a1_only_directed_contact",
    "a2_only_directed_contact",
    "any_social_contact"
  )

  purrr::map_dfr(event_cols, function(ev) {
    tmp <- summarise_event(frame_tbl[[ev]], fps, config$min_event_duration_frames)
    tmp |>
      mutate(
        file = unique(frame_tbl$file),
        event_type = ev,
        .before = 1
      )
  })
}

make_file_summary <- function(frame_tbl, event_summary_long, qc_tbl, config) {
  file_id <- unique(frame_tbl$file)

  a1_total <- event_summary_long |>
    filter(event_type == "a1_any_directed_contact") |>
    pull(duration_s)

  a2_total <- event_summary_long |>
    filter(event_type == "a2_any_directed_contact") |>
    pull(duration_s)

  any_contact <- event_summary_long |>
    filter(event_type == "any_social_contact") |>
    pull(duration_s)

  prox_close <- event_summary_long |>
    filter(event_type == "close_proximity") |>
    pull(duration_s)

  face_total <- event_summary_long |>
    filter(event_type == "face_investigation") |>
    pull(duration_s)

  body_total <- event_summary_long |>
    filter(event_type == "body_investigation") |>
    pull(duration_s)

  anogenital_total <- event_summary_long |>
    filter(event_type == "anogenital_investigation") |>
    pull(duration_s)

  mutual_total <- event_summary_long |>
    filter(event_type == "mutual_directed_contact") |>
    pull(duration_s)

  a1_only_total <- event_summary_long |>
    filter(event_type == "a1_only_directed_contact") |>
    pull(duration_s)

  a2_only_total <- event_summary_long |>
    filter(event_type == "a2_only_directed_contact") |>
    pull(duration_s)

  approach_total <- event_summary_long |>
    filter(event_type == "approach_event") |>
    pull(duration_s)

  retreat_total <- event_summary_long |>
    filter(event_type == "retreat_event") |>
    pull(duration_s)

  total_distance_moved <- sum(frame_tbl$speed1 / config$fps, na.rm = TRUE) +
    sum(frame_tbl$speed2 / config$fps, na.rm = TRUE)

  # Positive = more directed contact from animal 1 than animal 2.
  directional_bias <- safe_divide(a1_total - a2_total, a1_total + a2_total)

  tibble(
    file = file_id,
    n_frames = nrow(frame_tbl),
    duration_video_s = nrow(frame_tbl) / config$fps,
    duration_analyzed_s = sum(!is.na(frame_tbl$d_body_body)) / config$fps,
    any_social_contact_s = any_contact,
    close_proximity_s = prox_close,
    a1_directed_contact_s = a1_total,
    a2_directed_contact_s = a2_total,
    directional_bias_a1_minus_a2 = directional_bias,
    a1_minus_a2_directed_contact_s = a1_total - a2_total,
    a1_over_a2_directed_contact_ratio = safe_divide(a1_total, a2_total),
    a1_fraction_of_total_directed_contact = safe_divide(a1_total, a1_total + a2_total),
    face_investigation_s = face_total,
    body_investigation_s = body_total,
    anogenital_investigation_s = anogenital_total,
    mutual_directed_contact_s = mutual_total,
    a1_only_directed_contact_s = a1_only_total,
    a2_only_directed_contact_s = a2_only_total,
    approach_event_s = approach_total,
    retreat_event_s = retreat_total,
    total_distance_moved = total_distance_moved,
    directed_contact_per_distance_moved = safe_divide(a1_total + a2_total, total_distance_moved),
    close_proximity_per_distance_moved = safe_divide(prox_close, total_distance_moved),
    mean_body_distance = mean(frame_tbl$d_body_body, na.rm = TRUE),
    median_body_distance = median(frame_tbl$d_body_body, na.rm = TRUE),
    mean_speed_a1 = mean(frame_tbl$speed1, na.rm = TRUE),
    mean_speed_a2 = mean(frame_tbl$speed2, na.rm = TRUE),
    qc_max_missing_fraction = max(qc_tbl$frac_missing_original, na.rm = TRUE),
    qc_any_high_missing = any(qc_tbl$flag_high_missing, na.rm = TRUE),
    qc_total_impossible_jumps = sum(qc_tbl$impossible_jump_n, na.rm = TRUE),
    qc_fail = any(qc_tbl$qc_fail, na.rm = TRUE)
  )
}

plot_distance_events <- function(frame_tbl, file_id, output_dir, config) {
  p <- ggplot(frame_tbl, aes(x = time_s)) +
    geom_line(aes(y = d_body_body), linewidth = 0.3, alpha = 0.8) +
    geom_point(
      data = frame_tbl |> filter(any_social_contact),
      aes(y = d_body_body),
      size = 0.4,
      alpha = 0.8
    ) +
    labs(
      title = paste0(file_id, " — body-centre distance and social contact frames"),
      x = "Time (s)",
      y = paste0("Body-centre distance (", get_output_unit(config), ")")
    ) +
    theme_classic(base_size = 9)

  ggsave(
    filename = file.path(output_dir, "figures", paste0(file_id, "_distance_events.png")),
    plot = p,
    width = config$plot_width,
    height = config$plot_height,
    dpi = 300
  )
}

plot_trajectory <- function(Tracking, file_id, output_dir, config) {
  a1 <- get_xy(Tracking, "bodycentre_1") |>
    mutate(frame = row_number(), animal = "animal1")
  a2 <- get_xy(Tracking, "bodycentre_2") |>
    mutate(frame = row_number(), animal = "animal2")

  traj <- bind_rows(a1, a2)

  p <- ggplot(traj, aes(x = x, y = y, group = animal)) +
    geom_path(linewidth = 0.3, alpha = 0.7) +
    coord_equal() +
    labs(
      title = paste0(file_id, " — body-centre trajectories"),
      x = paste0("x (", get_output_unit(config), ")"),
      y = paste0("y (", get_output_unit(config), ")")
    ) +
    theme_classic(base_size = 9)

  ggsave(
    filename = file.path(output_dir, "figures", paste0(file_id, "_trajectories.png")),
    plot = p,
    width = config$plot_width,
    height = config$plot_height,
    dpi = 300
  )
}

analyze_socint_file <- function(input_file, config, bodyparts_all) {
  file_id <- tools::file_path_sans_ext(basename(input_file))
  message("Processing: ", file_id)

  summary_path <- file.path(config$output_dir, "tables", paste0(file_id, "_event_summary_long.csv"))
  if (file.exists(summary_path) && !isTRUE(config$overwrite_outputs)) {
    message("Skipping existing output for ", file_id)
    return(NULL)
  }

  Tracking <- ReadDLCDataFromCSV(file = input_file, fps = config$fps)
  validate_tracking(Tracking, bodyparts_all, config, file_id)

  calibration_res <- maybe_calibrate_tracking(Tracking, file_id, config, bodyparts_all)
  Tracking <- calibration_res$Tracking

  qc_res <- qc_tracking(Tracking, bodyparts_all, config, file_id)
  Tracking <- qc_res$Tracking
  qc_tbl <- qc_res$qc
  validate_qc(qc_tbl, config, file_id)

  frame_tbl <- make_event_table(Tracking, config, file_id)
  event_summary_long <- summarise_events_long(frame_tbl, config)
  file_summary <- make_file_summary(frame_tbl, event_summary_long, qc_tbl, config)
  proximity_bin_summary <- make_proximity_bin_summary(frame_tbl, config)
  overlap_summary <- make_overlap_summary(frame_tbl)
  interbout_summary <- make_interbout_summary(frame_tbl, config)

  file_summary <- file_summary |>
    left_join(proximity_bin_summary, by = "file") |>
    left_join(overlap_summary, by = "file") |>
    left_join(interbout_summary, by = "file")

  readr::write_csv(
    frame_tbl,
    file.path(config$output_dir, "frames", paste0(file_id, "_frame_events.csv"))
  )

  readr::write_csv(
    event_summary_long,
    file.path(config$output_dir, "tables", paste0(file_id, "_event_summary_long.csv"))
  )

  readr::write_csv(
    qc_tbl,
    file.path(config$output_dir, "qc", paste0(file_id, "_tracking_qc.csv"))
  )

  readr::write_csv(
    calibration_res$calibration,
    file.path(config$output_dir, "qc", paste0(file_id, "_arena_calibration.csv"))
  )

  if (nrow(calibration_res$arena_geom) > 0) {
    readr::write_csv(
      calibration_res$arena_geom,
      file.path(config$output_dir, "qc", paste0(file_id, "_arena_geometry.csv"))
    )
  }

  if (isTRUE(config$save_plots)) {
    plot_distance_events(frame_tbl, file_id, config$output_dir, config)
    plot_trajectory(Tracking, file_id, config$output_dir, config)
  }

  list(
    file_summary = file_summary,
    event_summary_long = event_summary_long,
    qc = qc_tbl,
    frame_tbl = frame_tbl,
    calibration = calibration_res$calibration,
    arena_geom = calibration_res$arena_geom
  )
}

# -------------------------------
# Additional summary functions
# -------------------------------

make_proximity_bin_summary <- function(frame_tbl, config) {
  bins <- config$proximity_bins
  labels <- config$proximity_bin_labels

  if (length(labels) != length(bins) - 1) {
    stop("proximity_bin_labels must have length length(proximity_bins) - 1")
  }

  frame_tbl |>
    mutate(
      proximity_bin = cut(
        d_body_body,
        breaks = bins,
        labels = labels,
        include.lowest = TRUE,
        right = FALSE
      )
    ) |>
    filter(!is.na(proximity_bin)) |>
    count(file, proximity_bin, name = "frames") |>
    group_by(file) |>
    mutate(
      total_valid_frames = sum(frames),
      time_s = frames / config$fps,
      percent_time = safe_divide(frames, total_valid_frames) * 100
    ) |>
    ungroup() |>
    select(file, proximity_bin, time_s, percent_time) |>
    tidyr::pivot_wider(
      names_from = proximity_bin,
      values_from = c(time_s, percent_time),
      names_glue = "proximity_{proximity_bin}_{.value}",
      values_fill = 0
    )
}

make_overlap_summary <- function(frame_tbl) {
  directed <- frame_tbl$a1_any_directed_contact | frame_tbl$a2_any_directed_contact
  close <- frame_tbl$close_proximity
  any_contact <- frame_tbl$any_social_contact
  mutual <- frame_tbl$mutual_directed_contact

  tibble(
    file = unique(frame_tbl$file),
    fraction_directed_contact_during_close_proximity = safe_divide(sum(directed & close, na.rm = TRUE), sum(close, na.rm = TRUE)),
    fraction_any_contact_during_close_proximity = safe_divide(sum(any_contact & close, na.rm = TRUE), sum(close, na.rm = TRUE)),
    fraction_mutual_contact_during_directed_contact = safe_divide(sum(mutual & directed, na.rm = TRUE), sum(directed, na.rm = TRUE))
  )
}

make_interbout_summary <- function(frame_tbl, config) {
  intervals <- interbout_intervals_s(frame_tbl$any_social_contact, config$fps, config$min_event_duration_frames)

  tibble(
    file = unique(frame_tbl$file),
    any_social_contact_bout_rate_per_min = summarise_event(frame_tbl$any_social_contact, config$fps, config$min_event_duration_frames)$bout_n / (nrow(frame_tbl) / config$fps / 60),
    any_social_contact_mean_interbout_interval_s = ifelse(length(intervals) == 0, NA_real_, mean(intervals)),
    any_social_contact_median_interbout_interval_s = ifelse(length(intervals) == 0, NA_real_, median(intervals))
  )
}

# -------------------------------
# 5) Run batch analysis
# -------------------------------

file_list <- list.files(
  path = config$input_dir,
  pattern = config$csv_pattern,
  recursive = isTRUE(config$recursive_input_search),
  full.names = TRUE
)

if (length(file_list) == 0) {
  stop("No CSV files found in input_dir: ", config$input_dir)
}

failures <- list()
results <- purrr::map(file_list, function(input_file) {
  tryCatch(
    analyze_socint_file(input_file, config = config, bodyparts_all = bodyparts_all),
    error = function(e) {
      warning("Failed to process ", basename(input_file), ": ", conditionMessage(e))
      failures[[length(failures) + 1]] <<- tibble(
        file = tools::file_path_sans_ext(basename(input_file)),
        input_file = input_file,
        error = conditionMessage(e)
      )
      NULL
    }
  )
})
results <- purrr::compact(results)

if (length(results) == 0) {
  stop("No files were successfully processed. Check warnings and input files.")
}

combined_file_summary <- purrr::map_dfr(results, "file_summary")
combined_event_summary_long <- purrr::map_dfr(results, "event_summary_long")
combined_qc <- purrr::map_dfr(results, "qc")
combined_calibration <- purrr::map_dfr(results, "calibration")
combined_arena_geometry <- purrr::map_dfr(results, "arena_geom")

# Wide version of event summary: useful for statistics.
combined_event_summary_wide <- combined_event_summary_long |>
  select(file, event_type, duration_s, percent_time, bout_n, latency_s, mean_bout_s, max_bout_s) |>
  pivot_wider(
    names_from = event_type,
    values_from = c(duration_s, percent_time, bout_n, latency_s, mean_bout_s, max_bout_s),
    names_glue = "{event_type}_{.value}"
  )

combined_output <- combined_file_summary |>
  left_join(combined_event_summary_wide, by = "file")

readr::write_csv(
  combined_file_summary,
  file.path(config$output_dir, "tables", "combined_file_summary.csv")
)

readr::write_csv(
  combined_event_summary_long,
  file.path(config$output_dir, "tables", "combined_event_summary_long.csv")
)

readr::write_csv(
  combined_event_summary_wide,
  file.path(config$output_dir, "tables", "combined_event_summary_wide.csv")
)

readr::write_csv(
  combined_qc,
  file.path(config$output_dir, "qc", "combined_tracking_qc.csv")
)

readr::write_csv(
  combined_calibration,
  file.path(config$output_dir, "qc", "combined_arena_calibration.csv")
)

if (nrow(combined_arena_geometry) > 0) {
  readr::write_csv(
    combined_arena_geometry,
    file.path(config$output_dir, "qc", "combined_arena_geometry.csv")
  )
}

readr::write_csv(
  combined_output,
  file.path(config$output_dir, "combined_output_event_based.csv")
)


if (length(failures) > 0) {
  failure_tbl <- bind_rows(failures)
  readr::write_csv(
    failure_tbl,
    file.path(config$output_dir, "tables", "failures.csv")
  )
  warning(nrow(failure_tbl), " file(s) failed. See tables/failures.csv.")
}

message("Done. Output written to: ", config$output_dir)

# -------------------------------
# 6) Notes for validation
# -------------------------------
#
# Recommended checks:
# 1. Open several *_distance_events.png plots and confirm event dots occur during real contact/proximity.
# 2. Open *_frame_events.csv for a few videos and inspect frames where any_social_contact == TRUE.
# 3. Check combined_tracking_qc.csv before trusting metrics.
# 4. Tune thresholds:
#      contact_dist
#      side_by_side_dist
#      close_proximity_dist
#      medium_proximity_dist
#      facing_angle_deg
#    after comparing against manually scored example videos.
# 5. If coordinates are calibrated to cm, rename unit = "cm" and adjust thresholds accordingly.
#
# ================================================================
