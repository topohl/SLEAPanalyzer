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

  # Set these paths.
  # Use forward slashes also on Windows to avoid escaping problems.
  working_dir = "C:/Users/topohl/Documents/GitHub/SLEAPanalyzer/SLEAPanalzyer",
  functions_file = "DLCAnalyzer_Functions_final.R",

  input_dir  = "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC",
  output_dir = "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC/output_event_based",

  # File pattern.
  csv_pattern = "\\.csv$",

  # Coordinate units.
  # If data are not calibrated, these are pixels.
  # If data are calibrated to cm, set thresholds in cm.
  unit = "px",

  # Interaction thresholds.
  contact_dist = 30,
  side_by_side_dist = 80,
  close_proximity_dist = 80,
  medium_proximity_dist = 160,

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

if (!dir.exists(config$working_dir)) {
  warning("working_dir does not exist: ", config$working_dir)
} else {
  setwd(config$working_dir)
}

if (file.exists(config$functions_file)) {
  source(config$functions_file)
} else {
  stop("Could not find functions file: ", config$functions_file,
       "\nExpected function ReadDLCDataFromCSV(file, fps).")
}

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
  ifelse(is.na(y) | y == 0, NA_real_, x / y)
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
  ifelse(is.na(idx), NA_real_, (idx - 1) / fps)
}

summarise_event <- function(event_vec, fps) {
  event_vec <- ifelse(is.na(event_vec), FALSE, event_vec)

  starts <- event_vec & !dplyr::lag(event_vec, default = FALSE)
  ends   <- event_vec & !dplyr::lead(event_vec, default = FALSE)

  bout_id <- cumsum(starts)
  bout_id[!event_vec] <- NA_integer_

  bout_df <- tibble(active = event_vec, bout_id = bout_id) |>
    filter(active) |>
    group_by(bout_id) |>
    summarise(frames = n(), .groups = "drop")

  tibble(
    duration_s = sum(event_vec) / fps,
    percent_time = 100 * mean(event_vec),
    bout_n = sum(starts),
    latency_s = first_true_time(event_vec, fps),
    mean_bout_s = ifelse(nrow(bout_df) == 0, NA_real_, mean(bout_df$frames) / fps),
    max_bout_s  = ifelse(nrow(bout_df) == 0, NA_real_, max(bout_df$frames) / fps),
    event_frames = sum(event_vec)
  )
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
        flag_missing_bodypart = TRUE
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
        flag_missing_bodypart = FALSE
      )
    }
  }

  list(
    Tracking = Tracking,
    qc = bind_rows(qc_rows)
  )
}

make_event_table <- function(Tracking, config, file_id) {
  fps <- config$fps
  n_frames <- length(Tracking$data$nose_1$x)

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

      a1_nose_to_nose2 = a1_nose_nose,
      a1_nose_to_body2 = a1_nose_body,
      a1_nose_to_tail2 = a1_nose_tail,
      a2_nose_to_nose1 = a2_nose_nose,
      a2_nose_to_body1 = a2_nose_body,
      a2_nose_to_tail1 = a2_nose_tail,

      side_by_side_parallel = side_by_side_parallel,
      side_by_side_antiparallel = side_by_side_antiparallel,
      close_proximity = close_proximity,
      medium_proximity = medium_proximity,
      mutually_facing = mutually_facing,
      a1_following_a2 = a1_following_a2,
      a2_following_a1 = a2_following_a1,
      a1_avoidance_from_a2 = a1_avoidance_from_a2,
      a2_avoidance_from_a1 = a2_avoidance_from_a1,

      # Combined directional investigation states.
      a1_any_directed_contact = a1_nose_nose | a1_nose_body | a1_nose_tail,
      a2_any_directed_contact = a2_nose_nose | a2_nose_body | a2_nose_tail,
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
    "side_by_side_parallel",
    "side_by_side_antiparallel",
    "close_proximity",
    "medium_proximity",
    "mutually_facing",
    "a1_following_a2",
    "a2_following_a1",
    "a1_avoidance_from_a2",
    "a2_avoidance_from_a1",
    "a1_any_directed_contact",
    "a2_any_directed_contact",
    "any_social_contact"
  )

  purrr::map_dfr(event_cols, function(ev) {
    tmp <- summarise_event(frame_tbl[[ev]], fps)
    tmp |>
      mutate(
        file = unique(frame_tbl$file),
        event_type = ev,
        .before = 1
      )
  })
}

make_file_summary <- function(frame_tbl, event_summary_long, qc_tbl) {
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

  # Positive = more directed contact from animal 1 than animal 2.
  directional_bias <- safe_divide(a1_total - a2_total, a1_total + a2_total)

  tibble(
    file = file_id,
    n_frames = nrow(frame_tbl),
    duration_video_s = max(frame_tbl$time_s, na.rm = TRUE),
    any_social_contact_s = any_contact,
    close_proximity_s = prox_close,
    a1_directed_contact_s = a1_total,
    a2_directed_contact_s = a2_total,
    directional_bias_a1_minus_a2 = directional_bias,
    mean_body_distance = mean(frame_tbl$d_body_body, na.rm = TRUE),
    median_body_distance = median(frame_tbl$d_body_body, na.rm = TRUE),
    mean_speed_a1 = mean(frame_tbl$speed1, na.rm = TRUE),
    mean_speed_a2 = mean(frame_tbl$speed2, na.rm = TRUE),
    qc_max_missing_fraction = max(qc_tbl$frac_missing_original, na.rm = TRUE),
    qc_any_high_missing = any(qc_tbl$flag_high_missing, na.rm = TRUE),
    qc_total_impossible_jumps = sum(qc_tbl$impossible_jump_n, na.rm = TRUE)
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
      y = paste0("Body-centre distance (", config$unit, ")")
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
      x = paste0("x (", config$unit, ")"),
      y = paste0("y (", config$unit, ")")
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

  Tracking <- ReadDLCDataFromCSV(file = input_file, fps = config$fps)

  missing_bodyparts <- setdiff(bodyparts_all, names(Tracking$data))
  if (length(missing_bodyparts) > 0) {
    warning("Missing bodyparts in ", file_id, ": ", paste(missing_bodyparts, collapse = ", "))
  }

  qc_res <- qc_tracking(Tracking, bodyparts_all, config, file_id)
  Tracking <- qc_res$Tracking
  qc_tbl <- qc_res$qc

  frame_tbl <- make_event_table(Tracking, config, file_id)
  event_summary_long <- summarise_events_long(frame_tbl, config)
  file_summary <- make_file_summary(frame_tbl, event_summary_long, qc_tbl)

  # Save per-file tables.
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

  if (isTRUE(config$save_plots)) {
    plot_distance_events(frame_tbl, file_id, config$output_dir, config)
    plot_trajectory(Tracking, file_id, config$output_dir, config)
  }

  list(
    file_summary = file_summary,
    event_summary_long = event_summary_long,
    qc = qc_tbl,
    frame_tbl = frame_tbl
  )
}

# -------------------------------
# 5) Run batch analysis
# -------------------------------

file_list <- list.files(
  path = config$input_dir,
  pattern = config$csv_pattern,
  full.names = TRUE
)

if (length(file_list) == 0) {
  stop("No CSV files found in input_dir: ", config$input_dir)
}

results <- purrr::map(
  file_list,
  ~ analyze_socint_file(.x, config = config, bodyparts_all = bodyparts_all)
)

combined_file_summary <- purrr::map_dfr(results, "file_summary")
combined_event_summary_long <- purrr::map_dfr(results, "event_summary_long")
combined_qc <- purrr::map_dfr(results, "qc")

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
  combined_output,
  file.path(config$output_dir, "combined_output_event_based.csv")
)

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
