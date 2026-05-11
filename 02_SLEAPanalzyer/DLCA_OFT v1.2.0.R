# ================================================================
# SLEAP/DLC Open Field Test Analyzer
# QC-first OFT analysis with frame-level, summary, and enhanced outputs
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
# Enhanced metrics now integrated directly into v1.2.0:
# - Habituation slopes across time bins
# - First-vs-last epoch summaries and deltas
# - Wall-threshold thigmotaxis metrics
# - Spatial occupancy entropy and arena coverage
# - Path tortuosity and angular movement metrics
# - Velocity-state segmentation
# - Composite center-exploration profile score across processed animals
# - Nature-oriented SVG summary plots
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
  min_valid_frames_fraction = 0.70,

  # Movement analysis.
  movement_cutoff_cm_s = 5,
  integration_period_frames = 5,
  immobility_min_duration_s = 1.0,

  # Time-course output.
  bin_size_s = 60,

  # Enhanced OFT metrics.
  epoch_duration_s = 300,
  wall_thresholds_cm = c(5, 10),
  occupancy_grid_n = 10,
  speed_low_cm_s = 2,
  speed_high_cm_s = 10,

  # Plots.
  save_plots = TRUE,
  save_enhanced_plots = TRUE,
  plot_width = 7,
  plot_height = 4,
  plot_width_single_col = 85 / 25.4,
  plot_height_single_col = 65 / 25.4,
  plot_width_double_col = 180 / 25.4,
  plot_height_double_col = 85 / 25.4
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

safe_divide <- function(x, y) {
  ifelse(is.na(y) | y == 0, NA_real_, x / y)
}

zscore <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  sx <- stats::sd(x, na.rm = TRUE)
  mx <- mean(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(NA_real_, length(x)))
  (x - mx) / sx
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

fit_slope <- function(data, y, x = "bin_s") {
  dat <- data %>%
    filter(!is.na(.data[[y]]), !is.na(.data[[x]]))

  if (nrow(dat) < 3 || dplyr::n_distinct(dat[[x]]) < 3) {
    return(tibble(
      metric = y,
      slope_per_min = NA_real_,
      intercept = NA_real_,
      r_squared = NA_real_,
      p_value = NA_real_
    ))
  }

  mod <- stats::lm(stats::as.formula(paste(y, "~", x)), data = dat)
  sm <- summary(mod)

  tibble(
    metric = y,
    slope_per_min = unname(stats::coef(mod)[[x]]) * 60,
    intercept = unname(stats::coef(mod)[["(Intercept)"]]),
    r_squared = sm$r.squared,
    p_value = sm$coefficients[x, "Pr(>|t|)"]
  )
}

make_habituation_slopes <- function(binned_tbl) {
  slope_metrics <- c(
    "distance_cm", "mean_speed_cm_s", "center_time_s", "periphery_time_s",
    "corner_time_s", "moving_time_s", "immobile_time_s", "mean_wall_distance_cm",
    "center_entries"
  )

  purrr::map_dfr(slope_metrics, fit_slope, data = binned_tbl) %>%
    pivot_wider(
      names_from = metric,
      values_from = c(slope_per_min, intercept, r_squared, p_value),
      names_glue = "{.value}_{metric}"
    )
}

make_epoch_metrics <- function(frame_tbl, config) {
  max_time <- max(frame_tbl$time_s, na.rm = TRUE)
  epoch_s <- min(config$epoch_duration_s, max_time / 2)

  epoch_tbl <- frame_tbl %>%
    mutate(
      epoch = case_when(
        time_s < epoch_s ~ "first_epoch",
        time_s >= max_time - epoch_s ~ "last_epoch",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(epoch)) %>%
    group_by(file, Batch, Code, epoch) %>%
    summarise(
      epoch_time_s = n() / config$fps,
      distance_cm = sum(speed_cm_s / config$fps, na.rm = TRUE),
      mean_speed_cm_s = mean(speed_cm_s, na.rm = TRUE),
      center_time_s = sum(in_center, na.rm = TRUE) / config$fps,
      center_time_percent = safe_divide(center_time_s, epoch_time_s) * 100,
      periphery_time_s = sum(in_periphery, na.rm = TRUE) / config$fps,
      corner_time_s = sum(in_corners, na.rm = TRUE) / config$fps,
      moving_time_s = sum(is_moving, na.rm = TRUE) / config$fps,
      immobile_time_s = sum(is_immobile, na.rm = TRUE) / config$fps,
      mean_wall_distance_cm = mean(wall_distance_cm, na.rm = TRUE),
      wall_5cm_time_s = sum(wall_distance_cm <= 5, na.rm = TRUE) / config$fps,
      wall_10cm_time_s = sum(wall_distance_cm <= 10, na.rm = TRUE) / config$fps,
      wall_5cm_percent = safe_divide(wall_5cm_time_s, epoch_time_s) * 100,
      wall_10cm_percent = safe_divide(wall_10cm_time_s, epoch_time_s) * 100,
      .groups = "drop"
    )

  delta_tbl <- epoch_tbl %>%
    select(file, Batch, Code, epoch, distance_cm, mean_speed_cm_s, center_time_s,
           center_time_percent, mean_wall_distance_cm, wall_5cm_percent, wall_10cm_percent,
           moving_time_s, immobile_time_s) %>%
    pivot_wider(
      names_from = epoch,
      values_from = c(distance_cm, mean_speed_cm_s, center_time_s, center_time_percent,
                      mean_wall_distance_cm, wall_5cm_percent, wall_10cm_percent,
                      moving_time_s, immobile_time_s)
    ) %>%
    mutate(
      delta_distance_cm = distance_cm_last_epoch - distance_cm_first_epoch,
      delta_mean_speed_cm_s = mean_speed_cm_s_last_epoch - mean_speed_cm_s_first_epoch,
      delta_center_time_s = center_time_s_last_epoch - center_time_s_first_epoch,
      delta_center_time_percent = center_time_percent_last_epoch - center_time_percent_first_epoch,
      delta_mean_wall_distance_cm = mean_wall_distance_cm_last_epoch - mean_wall_distance_cm_first_epoch,
      delta_wall_5cm_percent = wall_5cm_percent_last_epoch - wall_5cm_percent_first_epoch,
      delta_wall_10cm_percent = wall_10cm_percent_last_epoch - wall_10cm_percent_first_epoch,
      delta_moving_time_s = moving_time_s_last_epoch - moving_time_s_first_epoch,
      delta_immobile_time_s = immobile_time_s_last_epoch - immobile_time_s_first_epoch
    )

  list(epoch = epoch_tbl, delta = delta_tbl)
}

make_wall_threshold_metrics <- function(frame_tbl, config) {
  total_time_s <- nrow(frame_tbl) / config$fps

  purrr::map_dfr(config$wall_thresholds_cm, function(thr) {
    near_wall <- tidyr::replace_na(frame_tbl$wall_distance_cm <= thr, FALSE)
    tibble(
      threshold_cm = thr,
      wall_time_s = sum(near_wall, na.rm = TRUE) / config$fps,
      wall_time_percent = safe_divide(wall_time_s, total_time_s) * 100,
      wall_entries = entry_count(near_wall)
    )
  }) %>%
    pivot_wider(
      names_from = threshold_cm,
      values_from = c(wall_time_s, wall_time_percent, wall_entries),
      names_glue = "{.value}_{threshold_cm}cm"
    )
}

make_velocity_state_metrics <- function(frame_tbl, config) {
  frame_tbl %>%
    mutate(
      velocity_state = case_when(
        is.na(speed_cm_s) ~ NA_character_,
        speed_cm_s < config$speed_low_cm_s ~ "low_speed_or_rest",
        speed_cm_s < config$speed_high_cm_s ~ "moderate_speed",
        TRUE ~ "high_speed"
      )
    ) %>%
    filter(!is.na(velocity_state)) %>%
    group_by(velocity_state) %>%
    summarise(
      state_time_s = n() / config$fps,
      state_mean_speed_cm_s = mean(speed_cm_s, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(state_time_percent = safe_divide(state_time_s, sum(state_time_s, na.rm = TRUE)) * 100) %>%
    pivot_wider(
      names_from = velocity_state,
      values_from = c(state_time_s, state_time_percent, state_mean_speed_cm_s),
      names_glue = "{.value}_{velocity_state}"
    )
}

make_path_geometry_metrics <- function(frame_tbl, config) {
  dat <- frame_tbl %>%
    arrange(frame) %>%
    mutate(
      dx = x - dplyr::lag(x),
      dy = y - dplyr::lag(y),
      step_cm = sqrt(dx^2 + dy^2),
      heading_rad = atan2(dy, dx),
      turn_angle_rad = atan2(sin(heading_rad - dplyr::lag(heading_rad)), cos(heading_rad - dplyr::lag(heading_rad))),
      turn_angle_deg = abs(turn_angle_rad) * 180 / pi,
      angular_velocity_deg_s = turn_angle_deg * config$fps
    )

  valid_xy <- dat %>% filter(!is.na(x), !is.na(y))
  net_displacement <- if (nrow(valid_xy) >= 2) {
    sqrt((dplyr::last(valid_xy$x) - dplyr::first(valid_xy$x))^2 +
           (dplyr::last(valid_xy$y) - dplyr::first(valid_xy$y))^2)
  } else {
    NA_real_
  }

  total_distance <- sum(dat$step_cm, na.rm = TRUE)

  tibble(
    valid_frame_fraction = mean(!is.na(dat$x) & !is.na(dat$y)),
    total_distance_from_frames_cm = total_distance,
    net_displacement_cm = net_displacement,
    tortuosity_distance_over_displacement = safe_divide(total_distance, net_displacement),
    mean_turn_angle_deg = mean(dat$turn_angle_deg, na.rm = TRUE),
    median_turn_angle_deg = median(dat$turn_angle_deg, na.rm = TRUE),
    mean_angular_velocity_deg_s = mean(dat$angular_velocity_deg_s, na.rm = TRUE),
    high_turn_fraction_90deg = mean(dat$turn_angle_deg >= 90, na.rm = TRUE)
  )
}

make_occupancy_entropy_metrics <- function(frame_tbl, config) {
  grid_n <- config$occupancy_grid_n
  valid <- frame_tbl %>% filter(!is.na(x), !is.na(y))

  if (nrow(valid) == 0) {
    return(tibble(
      occupancy_grid_n = grid_n,
      occupied_bins_n = NA_integer_,
      total_bins_n = grid_n * grid_n,
      arena_coverage_percent = NA_real_,
      occupancy_entropy = NA_real_,
      occupancy_entropy_normalized = NA_real_
    ))
  }

  x_range <- range(valid$x, na.rm = TRUE)
  y_range <- range(valid$y, na.rm = TRUE)

  valid %>%
    mutate(
      x_scaled = ifelse(diff(x_range) == 0, 0, (x - x_range[1]) / diff(x_range)),
      y_scaled = ifelse(diff(y_range) == 0, 0, (y - y_range[1]) / diff(y_range)),
      x_bin = pmin(pmax(floor(x_scaled * grid_n), 0), grid_n - 1),
      y_bin = pmin(pmax(floor(y_scaled * grid_n), 0), grid_n - 1),
      spatial_bin = paste(x_bin, y_bin, sep = "_")
    ) %>%
    count(spatial_bin, name = "n_frames") %>%
    summarise(
      occupancy_grid_n = grid_n,
      occupied_bins_n = n(),
      total_bins_n = grid_n * grid_n,
      arena_coverage_percent = occupied_bins_n / total_bins_n * 100,
      occupancy_entropy = {
        p <- n_frames / sum(n_frames)
        -sum(p * log2(p), na.rm = TRUE)
      },
      occupancy_entropy_normalized = occupancy_entropy / log2(total_bins_n)
    )
}

make_enhanced_row <- function(summary_tbl, frame_tbl, binned_tbl, config) {
  epoch_metrics <- make_epoch_metrics(frame_tbl, config)

  bind_cols(
    summary_tbl %>% select(file, Batch, Code),
    make_habituation_slopes(binned_tbl),
    epoch_metrics$delta %>% select(-file, -Batch, -Code),
    make_wall_threshold_metrics(frame_tbl, config),
    make_velocity_state_metrics(frame_tbl, config),
    make_path_geometry_metrics(frame_tbl, config),
    make_occupancy_entropy_metrics(frame_tbl, config)
  )
}

add_center_exploration_score <- function(summary_tbl) {
  summary_tbl %>%
    mutate(
      z_center_time_percent = zscore(center_time_percent),
      z_center_entries = zscore(center_entries),
      z_center_latency_inverse = -zscore(center_latency_s),
      z_mean_wall_distance = zscore(mean_wall_distance_cm),
      oft_center_exploration_score = rowMeans(
        cbind(
          z_center_time_percent,
          z_center_entries,
          z_center_latency_inverse,
          z_mean_wall_distance
        ),
        na.rm = TRUE
      ),
      oft_center_exploration_score = ifelse(is.nan(oft_center_exploration_score), NA_real_, oft_center_exploration_score)
    )
}

make_theme_nature <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.text = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0)
    )
}

save_enhanced_plots <- function(summary_tbl, bins_tbl, output_dir, config) {
  if (!isTRUE(config$save_enhanced_plots)) return(invisible(NULL))
  if (nrow(summary_tbl) == 0) return(invisible(NULL))

  plot_dir <- file.path(output_dir, "enhanced_plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  if ("oft_center_exploration_score" %in% names(summary_tbl)) {
    p1 <- summary_tbl %>%
      ggplot(aes(x = reorder(file, oft_center_exploration_score), y = oft_center_exploration_score)) +
      geom_col(width = 0.75) +
      coord_flip() +
      labs(
        title = "OFT center-exploration profile score",
        x = NULL,
        y = "Composite z-score"
      ) +
      make_theme_nature()

    ggsave(
      file.path(plot_dir, "enhanced_oft_center_exploration_score.svg"),
      p1,
      width = config$plot_width_double_col,
      height = config$plot_height_double_col
    )
  }

  p2 <- summary_tbl %>%
    ggplot(aes(x = mean_wall_distance_cm, y = center_time_percent)) +
    geom_point(size = 1.8, alpha = 0.85) +
    labs(
      title = "Center occupancy versus wall distance",
      x = "Mean wall distance (cm)",
      y = "Center time (%)"
    ) +
    make_theme_nature()

  ggsave(
    file.path(plot_dir, "enhanced_center_vs_wall_distance.svg"),
    p2,
    width = config$plot_width_single_col,
    height = config$plot_height_single_col
  )

  if (nrow(bins_tbl) > 0) {
    p3 <- bins_tbl %>%
      ggplot(aes(x = bin_s / 60, y = distance_cm, group = file)) +
      geom_line(alpha = 0.35, linewidth = 0.3) +
      stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 0.8) +
      labs(
        title = "OFT locomotor habituation trajectory",
        x = "Time (min)",
        y = "Distance per bin (cm)"
      ) +
      make_theme_nature()

    ggsave(
      file.path(plot_dir, "enhanced_distance_habituation_trajectory.svg"),
      p3,
      width = config$plot_width_single_col,
      height = config$plot_height_single_col
    )

    p4 <- bins_tbl %>%
      ggplot(aes(x = bin_s / 60, y = center_time_s, group = file)) +
      geom_line(alpha = 0.35, linewidth = 0.3) +
      stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 0.8) +
      labs(
        title = "OFT center exploration over time",
        x = "Time (min)",
        y = "Center time per bin (s)"
      ) +
      make_theme_nature()

    ggsave(
      file.path(plot_dir, "enhanced_center_time_trajectory.svg"),
      p4,
      width = config$plot_width_single_col,
      height = config$plot_height_single_col
    )
  }

  invisible(NULL)
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
  enhanced_tbl <- make_enhanced_row(summary_tbl, frame_tbl, binned_tbl, config)
  epoch_tbl <- make_epoch_metrics(frame_tbl, config)$epoch

  readr::write_csv(summary_tbl, file.path(dirs$tables, paste0(file_id, "_summary.csv")))
  readr::write_csv(frame_tbl, file.path(dirs$frames, paste0(file_id, "_frames.csv")))
  readr::write_csv(binned_tbl, file.path(dirs$bins, paste0(file_id, "_bins.csv")))
  readr::write_csv(qc_tbl, file.path(dirs$qc, paste0(file_id, "_qc.csv")))
  readr::write_csv(enhanced_tbl, file.path(dirs$enhanced, paste0(file_id, "_enhanced_metrics.csv")))
  readr::write_csv(epoch_tbl, file.path(dirs$enhanced, paste0(file_id, "_first_last_epoch_metrics.csv")))

  save_optional_plots(tracking, file_id, dirs$plots, config)

  list(
    summary = summary_tbl,
    frames = frame_tbl,
    bins = binned_tbl,
    qc = qc_tbl,
    enhanced = enhanced_tbl,
    epochs = epoch_tbl
  )
}

# -------------------------------
# 4) Batch processing
# -------------------------------

all_summaries <- list()
all_bins <- list()
all_qc <- list()
all_enhanced <- list()
all_epochs <- list()
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
  dirs$enhanced <- file.path(dirs$output, "enhanced_metrics")

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
      all_enhanced[[length(all_enhanced) + 1]] <- result$enhanced
      all_epochs[[length(all_epochs) + 1]] <- result$epochs
    }
  }

  batch_summaries <- bind_rows(all_summaries) %>% filter(Batch == batch)
  batch_bins <- bind_rows(all_bins) %>% filter(Batch == batch)
  batch_qc <- bind_rows(all_qc) %>% filter(Batch == batch)
  batch_enhanced <- bind_rows(all_enhanced) %>% filter(Batch == batch)
  batch_epochs <- bind_rows(all_epochs) %>% filter(Batch == batch)

  if (nrow(batch_summaries) > 0) {
    batch_summaries_scored <- add_center_exploration_score(batch_summaries)
    batch_enhanced_scored <- batch_enhanced %>%
      left_join(
        batch_summaries_scored %>% select(file, Batch, Code, starts_with("z_"), oft_center_exploration_score),
        by = c("file", "Batch", "Code")
      ) %>%
      mutate(
        enhanced_qc_flag = batch_summaries$qc_flag[match(file, batch_summaries$file)] |
          valid_frame_fraction < config$min_valid_frames_fraction |
          is.na(occupancy_entropy_normalized)
      )

    readr::write_csv(batch_summaries_scored, file.path(dirs$output, "combined_summary.csv"))
    readr::write_csv(batch_bins, file.path(dirs$output, "combined_bins.csv"))
    readr::write_csv(batch_qc, file.path(dirs$output, "combined_qc.csv"))
    readr::write_csv(batch_enhanced_scored, file.path(dirs$output, "combined_enhanced_metrics.csv"))
    readr::write_csv(batch_epochs, file.path(dirs$output, "combined_first_last_epoch_metrics.csv"))

    save_enhanced_plots(batch_summaries_scored, batch_bins, dirs$output, config)
  }

  message("Processing complete for batch ", batch)
}

if (length(all_summaries) > 0) {
  global_output <- file.path(config$behavior_root, "OFT_SLEAP_output_v1.2.0")
  dir.create(global_output, recursive = TRUE, showWarnings = FALSE)

  all_summaries_tbl <- add_center_exploration_score(bind_rows(all_summaries))
  all_bins_tbl <- bind_rows(all_bins)
  all_qc_tbl <- bind_rows(all_qc)
  all_enhanced_tbl <- bind_rows(all_enhanced) %>%
    left_join(
      all_summaries_tbl %>% select(file, Batch, Code, starts_with("z_"), oft_center_exploration_score, qc_flag),
      by = c("file", "Batch", "Code")
    ) %>%
    mutate(
      enhanced_qc_flag = qc_flag |
        valid_frame_fraction < config$min_valid_frames_fraction |
        is.na(occupancy_entropy_normalized)
    )
  all_epochs_tbl <- bind_rows(all_epochs)

  readr::write_csv(all_summaries_tbl, file.path(global_output, "all_batches_summary.csv"))
  readr::write_csv(all_bins_tbl, file.path(global_output, "all_batches_bins.csv"))
  readr::write_csv(all_qc_tbl, file.path(global_output, "all_batches_qc.csv"))
  readr::write_csv(all_enhanced_tbl, file.path(global_output, "all_batches_enhanced_metrics.csv"))
  readr::write_csv(all_epochs_tbl, file.path(global_output, "all_batches_first_last_epoch_metrics.csv"))

  save_enhanced_plots(all_summaries_tbl, all_bins_tbl, global_output, config)
}

if (length(failures) > 0) {
  failure_tbl <- bind_rows(failures)
  failure_output <- file.path(config$behavior_root, "OFT_SLEAP_output_v1.2.0")
  dir.create(failure_output, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(failure_tbl, file.path(failure_output, "failures.csv"))
  warning(nrow(failure_tbl), " file(s) failed. See failures.csv.")
}

message("OFT analysis complete.")
