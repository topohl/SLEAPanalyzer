# Report-only quality summaries. These helpers never exclude or mutate frames.

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
      stringsAsFactors = FALSE
    )
  }))
}

tracking_qc_report <- function(tracking, landmarks = names(tracking$data),
                               likelihood_cutoff = NULL) {
  list(
    coordinates = coordinate_qc(tracking, landmarks),
    likelihood = likelihood_qc(tracking, landmarks, likelihood_cutoff),
    fps = get_tracking_fps(tracking),
    frames = length(get_tracking_frames(tracking)),
    duration_seconds = get_tracking_duration(tracking),
    coordinate_unit = get_tracking_unit(tracking)
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
