# Canonical tracking quality control.
#
# QC here is report-only: these helpers never silently exclude or mutate
# frames. They produce the numbers an analyst needs in order to decide whether
# a recording is usable, and every assay emits them alongside its results so a
# behavioral summary can always be read together with the quality of the
# tracking it rests on.

# ---------------------------------------------------------------------------
# Per-landmark summaries
# ---------------------------------------------------------------------------

coordinate_qc <- function(tracking, landmarks = names(tracking$data)) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  do.call(rbind, lapply(landmarks, function(point) {
    coordinates <- get_point_coordinates(tracking, point)
    valid <- is.finite(coordinates$x) & is.finite(coordinates$y)
    data.frame(
      landmark = point,
      frames = length(valid),
      valid_frames = sum(valid),
      missing_or_invalid_frames = sum(!valid),
      valid_fraction = mean(valid),
      stringsAsFactors = FALSE
    )
  }))
}

likelihood_qc <- function(tracking, landmarks = names(tracking$data), cutoff = NULL) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  if (!is.null(cutoff)) validate_scalar_number(cutoff, "cutoff")
  do.call(rbind, lapply(landmarks, function(point) {
    values <- tracking$data[[point]]$likelihood
    if (is.null(values)) values <- rep(NA_real_, length(tracking$frames))
    finite <- values[is.finite(values)]
    data.frame(
      landmark = point,
      observations = length(values),
      finite_observations = length(finite),
      mean_likelihood = if (length(finite) == 0) NA_real_ else mean(finite),
      min_likelihood = if (length(finite) == 0) NA_real_ else min(finite),
      below_cutoff = if (is.null(cutoff)) NA_integer_ else sum(values < cutoff, na.rm = TRUE),
      low_confidence_fraction = if (is.null(cutoff)) NA_real_ else
        mean(is.na(values) | values < cutoff),
      stringsAsFactors = FALSE
    )
  }))
}

#' Frames whose implied speed exceeds what the animal can achieve.
#'
#' The threshold is a speed in coordinate units per second, so it must be set
#' against calibrated data and against the species. A displacement flagged
#' here is usually an identity swap or a spurious detection, not locomotion.
#'
#' @param max_speed maximum plausible speed in coordinate units per second
displacement_qc <- function(tracking, landmarks = names(tracking$data),
                            max_speed = NULL) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  fps <- get_tracking_fps(tracking)
  if (!is.null(max_speed)) validate_scalar_number(max_speed, "max_speed", positive = TRUE)

  do.call(rbind, lapply(landmarks, function(point) {
    coordinates <- get_point_coordinates(tracking, point)
    dx <- c(NA_real_, diff(coordinates$x))
    dy <- c(NA_real_, diff(coordinates$y))
    speed <- sqrt(dx^2 + dy^2) * fps
    measurable <- is.finite(speed)
    data.frame(
      landmark = point,
      measurable_intervals = sum(measurable),
      max_speed = if (any(measurable)) max(speed[measurable]) else NA_real_,
      median_speed = if (any(measurable)) stats::median(speed[measurable]) else NA_real_,
      implausible_intervals = if (is.null(max_speed)) NA_integer_ else
        sum(measurable & speed > max_speed),
      implausible_fraction = if (is.null(max_speed) || !any(measurable)) NA_real_ else
        sum(measurable & speed > max_speed) / sum(measurable),
      stringsAsFactors = FALSE
    )
  }))
}

#' Distance between two landmarks that should be rigidly separated.
#'
#' A body length that varies far more than the animal can deform indicates a
#' mislabelled landmark or, in multi-animal recordings, an identity swap.
#'
#' @param front,back landmarks defining the body axis
#' @param tolerance fraction of the median length treated as acceptable
skeleton_qc <- function(tracking, front = "nose", back = "tailbase",
                        tolerance = 0.5) {
  validate_tracking_data(tracking, required_landmarks = c(front, back))
  validate_scalar_number(tolerance, "tolerance", positive = TRUE)
  length_series <- tracking_point_distance(tracking, front, back)
  measurable <- is.finite(length_series)
  if (!any(measurable)) {
    return(data.frame(
      front = front, back = back, measurable_frames = 0L,
      median_length = NA_real_, abnormal_frames = NA_integer_,
      abnormal_fraction = NA_real_, stringsAsFactors = FALSE
    ))
  }
  median_length <- stats::median(length_series[measurable])
  deviation <- abs(length_series - median_length) / median_length
  abnormal <- measurable & deviation > tolerance
  data.frame(
    front = front,
    back = back,
    measurable_frames = sum(measurable),
    median_length = median_length,
    abnormal_frames = sum(abnormal),
    abnormal_fraction = sum(abnormal) / sum(measurable),
    stringsAsFactors = FALSE
  )
}

#' Suspected identity swaps between two animals.
#'
#' A swap shows up as both animals appearing to jump simultaneously, with each
#' landing close to where the other was on the previous frame. Testing the
#' cross distance as well as the jump distinguishes a swap from two animals
#' that genuinely moved fast at the same moment.
#'
#' @param animal_a,animal_b the body-centre landmark of each animal
#' @param jump_threshold displacement, in coordinate units, treated as a jump
identity_swap_qc <- function(tracking, animal_a, animal_b, jump_threshold) {
  validate_tracking_data(tracking, required_landmarks = c(animal_a, animal_b))
  validate_scalar_number(jump_threshold, "jump_threshold", positive = TRUE)
  a <- get_point_coordinates(tracking, animal_a)
  b <- get_point_coordinates(tracking, animal_b)
  n <- nrow(a)
  if (n < 2) {
    return(data.frame(
      animal_a = animal_a, animal_b = animal_b, comparable_frames = 0L,
      suspected_swaps = 0L, suspected_swap_fraction = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  previous <- seq_len(n - 1L)
  current <- previous + 1L
  a_jump <- euclidean_distance(a$x[current], a$y[current], a$x[previous], a$y[previous])
  b_jump <- euclidean_distance(b$x[current], b$y[current], b$x[previous], b$y[previous])
  # Distance from each animal's new position to the other animal's old one.
  a_to_b_previous <- euclidean_distance(a$x[current], a$y[current], b$x[previous], b$y[previous])
  b_to_a_previous <- euclidean_distance(b$x[current], b$y[current], a$x[previous], a$y[previous])

  comparable <- is.finite(a_jump) & is.finite(b_jump) &
    is.finite(a_to_b_previous) & is.finite(b_to_a_previous)
  swap <- comparable &
    a_jump > jump_threshold & b_jump > jump_threshold &
    a_to_b_previous < a_jump & b_to_a_previous < b_jump

  data.frame(
    animal_a = animal_a,
    animal_b = animal_b,
    comparable_frames = sum(comparable),
    suspected_swaps = sum(swap),
    suspected_swap_fraction = if (any(comparable)) sum(swap) / sum(comparable) else NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Landmarks whose coordinate series disagree on length.
frame_consistency_qc <- function(tracking) {
  validate_tracking_data(tracking)
  n_frames <- length(get_tracking_frames(tracking))
  lengths <- vapply(tracking$data, function(point) nrow(point), integer(1))
  data.frame(
    frames = n_frames,
    landmarks = length(lengths),
    consistent = all(lengths == n_frames),
    min_landmark_frames = min(lengths),
    max_landmark_frames = max(lengths),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Combined report
# ---------------------------------------------------------------------------

#' One canonical QC report for a tracking file.
#'
#' `valid_time_s` is the analyzed time implied by the landmarks a measurement
#' depends on, and is the denominator every behavioral percentage should use.
#'
#' @param tracking a TrackingData object
#' @param landmarks landmarks to summarise
#' @param required_landmarks landmarks a measurement depends on; defaults to
#'   `landmarks`. Determines `valid_time_s`.
#' @param likelihood_cutoff minimum acceptable tracking confidence
#' @param max_speed maximum plausible speed, in coordinate units per second
#' @param skeleton a length-two character vector of body-axis landmarks
tracking_qc_report <- function(tracking, landmarks = names(tracking$data),
                               required_landmarks = landmarks,
                               likelihood_cutoff = NULL,
                               max_speed = NULL,
                               skeleton = NULL) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  fps <- get_tracking_fps(tracking)
  valid <- landmark_validity(
    tracking, required_landmarks,
    likelihood_cutoff = likelihood_cutoff,
    interpolated_is_valid = TRUE
  )
  observed_only <- landmark_validity(
    tracking, required_landmarks,
    likelihood_cutoff = likelihood_cutoff,
    interpolated_is_valid = FALSE
  )
  gaps <- rle(!valid)
  longest_gap <- if (any(gaps$values)) max(gaps$lengths[gaps$values]) else 0L

  list(
    file = tracking$filename,
    coordinates = coordinate_qc(tracking, landmarks),
    likelihood = likelihood_qc(tracking, landmarks, likelihood_cutoff),
    displacement = displacement_qc(tracking, landmarks, max_speed),
    provenance = if (any(vapply(tracking$data, function(p) !is.null(p$status), logical(1)))) {
      interpolation_report(tracking, landmarks)
    } else NULL,
    skeleton = if (is.null(skeleton)) NULL else skeleton_qc(tracking, skeleton[1], skeleton[2]),
    arena = if (is.null(tracking$arena_calibration)) NULL else
      arena_violation_report(tracking, landmarks),
    frames_consistent = frame_consistency_qc(tracking),
    fps = fps,
    frames = length(get_tracking_frames(tracking)),
    duration_seconds = get_tracking_duration(tracking),
    valid_frames = sum(valid),
    valid_time_s = sum(valid) / fps,
    valid_fraction = mean(valid),
    observed_time_s = sum(observed_only) / fps,
    interpolated_time_s = (sum(valid) - sum(observed_only)) / fps,
    longest_invalid_gap_frames = as.integer(longest_gap),
    longest_invalid_gap_s = longest_gap / fps,
    coordinate_unit = get_tracking_unit(tracking)
  )
}

#' Flatten a QC report to one row for a batch summary table.
qc_summary_row <- function(report) {
  data.frame(
    file = if (is.null(report$file)) NA_character_ else report$file,
    fps = report$fps,
    frames = report$frames,
    duration_s = report$duration_seconds,
    valid_time_s = report$valid_time_s,
    valid_fraction = report$valid_fraction,
    observed_time_s = report$observed_time_s,
    interpolated_time_s = report$interpolated_time_s,
    longest_invalid_gap_s = report$longest_invalid_gap_s,
    min_landmark_valid_fraction = min(report$coordinates$valid_fraction),
    coordinate_unit = report$coordinate_unit,
    frames_consistent = report$frames_consistent$consistent,
    stringsAsFactors = FALSE
  )
}

#' Apply QC thresholds and report which ones a recording fails.
#'
#' Returns the decision and the reasons; it does not drop anything. Excluding a
#' recording stays an explicit act by the analyst.
qc_flags <- function(report,
                     min_valid_fraction = 0.8,
                     max_longest_gap_s = 5,
                     max_interpolated_fraction = 0.2) {
  validate_scalar_number(min_valid_fraction, "min_valid_fraction", positive = TRUE, allow_zero = TRUE)
  validate_scalar_number(max_longest_gap_s, "max_longest_gap_s", positive = TRUE, allow_zero = TRUE)
  validate_scalar_number(max_interpolated_fraction, "max_interpolated_fraction", positive = TRUE, allow_zero = TRUE)

  interpolated_fraction <- if (report$duration_seconds > 0) {
    report$interpolated_time_s / report$duration_seconds
  } else NA_real_

  reasons <- character()
  if (isTRUE(report$valid_fraction < min_valid_fraction)) {
    reasons <- c(reasons, sprintf(
      "valid fraction %.3f below %.3f", report$valid_fraction, min_valid_fraction
    ))
  }
  if (isTRUE(report$longest_invalid_gap_s > max_longest_gap_s)) {
    reasons <- c(reasons, sprintf(
      "longest gap %.2f s above %.2f s", report$longest_invalid_gap_s, max_longest_gap_s
    ))
  }
  if (isTRUE(interpolated_fraction > max_interpolated_fraction)) {
    reasons <- c(reasons, sprintf(
      "interpolated fraction %.3f above %.3f", interpolated_fraction, max_interpolated_fraction
    ))
  }
  if (isFALSE(report$frames_consistent$consistent)) {
    reasons <- c(reasons, "landmark frame counts disagree")
  }

  list(
    pass = length(reasons) == 0,
    reasons = reasons,
    interpolated_fraction = interpolated_fraction
  )
}

geometry_qc_report <- function(corners) {
  tryCatch(
    rectangular_arena_geometry(corners),
    error = function(error) list(valid = FALSE, error = conditionMessage(error))
  )
}

calibration_qc_report <- function(corners, metric_width, metric_height) {
  tryCatch(
    c(
      list(valid = TRUE),
      rectangular_calibration_diagnostics(corners, metric_width, metric_height)
    ),
    error = function(error) list(valid = FALSE, error = conditionMessage(error))
  )
}
