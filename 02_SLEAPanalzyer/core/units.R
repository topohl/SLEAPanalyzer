# Explicit coordinate, threshold, time, and speed unit helpers.

canonical_coordinate_unit <- function(unit) {
  if (length(unit) != 1 || !is.character(unit) || is.na(unit) || !nzchar(unit)) {
    stop("coordinate unit must be one non-empty character value")
  }
  normalized <- tolower(trimws(unit))
  if (normalized %in% c("px", "pixel", "pixels")) return("px")
  if (normalized %in% c("cm", "centimeter", "centimeters", "centimetre", "centimetres")) return("cm")
  stop("Unsupported coordinate unit: ", unit, ". Expected pixels/px or cm.")
}

validate_coordinate_unit <- function(unit) {
  canonical_coordinate_unit(unit)
  unit
}

validate_threshold_unit <- function(threshold_unit, coordinate_unit) {
  threshold <- canonical_coordinate_unit(threshold_unit)
  coordinate <- canonical_coordinate_unit(coordinate_unit)
  if (!identical(threshold, coordinate)) {
    stop(
      "Threshold unit ('", threshold_unit,
      "') does not match coordinate unit ('", coordinate_unit, "')"
    )
  }
  threshold_unit
}

frames_to_seconds <- function(frames, fps) {
  validate_fps(fps)
  validate_numeric_vector(frames, "frames", finite = TRUE)
  frames / fps
}

#' Convert a duration in seconds to a whole number of frames.
#'
#' @param seconds a non-negative duration
#' @param fps frames per second
#' @param round_fn rounding applied to the frame count; use ceiling for
#'   minimum durations and floor for maximum tolerated gaps
seconds_to_frames <- function(seconds, fps, round_fn = ceiling) {
  validate_fps(fps)
  validate_scalar_number(seconds, "seconds", positive = TRUE, allow_zero = TRUE)
  if (!is.function(round_fn)) stop("round_fn must be a function")
  as.integer(round_fn(seconds * fps))
}

distance_to_speed <- function(distance, fps, interval_frames = 1,
                              coordinate_unit = NULL,
                              threshold_unit = coordinate_unit) {
  validate_fps(fps)
  validate_numeric_vector(distance, "distance")
  validate_scalar_number(interval_frames, "interval_frames", positive = TRUE)
  if (!is.null(coordinate_unit)) {
    validate_coordinate_unit(coordinate_unit)
    validate_threshold_unit(threshold_unit, coordinate_unit)
  }
  distance * fps / interval_frames
}

speed_unit <- function(coordinate_unit) {
  canonical <- canonical_coordinate_unit(coordinate_unit)
  paste0(canonical, "/s")
}
