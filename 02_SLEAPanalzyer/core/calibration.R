# Pure calibration calculations. Applying a scale to TrackingData remains explicit.

validate_calibration_scale <- function(scale) {
  validate_scalar_number(scale, "calibration scale", positive = TRUE)
  invisible(scale)
}

calibration_scale_from_distance <- function(pixel_distance, metric_distance) {
  validate_scalar_number(pixel_distance, "pixel_distance", positive = TRUE)
  validate_scalar_number(metric_distance, "metric_distance", positive = TRUE)
  metric_distance / pixel_distance
}

calibration_scale_from_area <- function(pixel_area, metric_area) {
  validate_scalar_number(pixel_area, "pixel_area", positive = TRUE)
  validate_scalar_number(metric_area, "metric_area", positive = TRUE)
  sqrt(metric_area / pixel_area)
}

rectangular_calibration_diagnostics <- function(corners, metric_width, metric_height) {
  geometry <- rectangular_arena_geometry(corners)
  if (!geometry$valid) stop("Arena geometry is degenerate and cannot be calibrated")
  validate_scalar_number(metric_width, "metric_width", positive = TRUE)
  validate_scalar_number(metric_height, "metric_height", positive = TRUE)
  horizontal_px <- mean(geometry$side_lengths[c(1, 3)])
  vertical_px <- mean(geometry$side_lengths[c(2, 4)])
  x_scale <- metric_width / horizontal_px
  y_scale <- metric_height / vertical_px
  list(
    x_scale = x_scale,
    y_scale = y_scale,
    mean_scale = mean(c(x_scale, y_scale)),
    anisotropy_ratio = max(x_scale, y_scale) / min(x_scale, y_scale),
    geometry = geometry
  )
}
