# Bounded gap interpolation with an explicit per-frame provenance record.
#
# The rule this module enforces: a coordinate is either observed, or it is a
# short interpolated bridge that is labelled as such, or it stays missing.
# Nothing else is acceptable, because a downstream measurement cannot tell a
# fabricated coordinate from a real one unless the pipeline records which is
# which.

FRAME_STATUS_LEVELS <- c("observed", "interpolated", "invalid")

#' Linearly interpolate short interior gaps only.
#'
#' Leading and trailing gaps are never filled: there is no observation on one
#' side to interpolate from, and forward/backward filling them invents a
#' stationary animal at the first or last known position. Interior gaps longer
#' than `max_gap_frames` are left missing rather than bridged across a period
#' where the animal may have moved anywhere.
#'
#' @param x a numeric vector
#' @param valid optional validity mask; combined with `is.finite(x)`
#' @param max_gap_frames longest interior gap that may be interpolated
#' @return a list with the filled `values` and a per-frame `status` vector
bounded_linear_interpolation <- function(x, valid = NULL, max_gap_frames = 0L) {
  validate_numeric_vector(x, "x")
  n <- length(x)
  valid <- normalize_validity_mask(valid, n) & is.finite(x)
  max_gap_frames <- validate_nonnegative_integer(max_gap_frames, "max_gap_frames")

  status <- ifelse(valid, "observed", "invalid")
  values <- x
  values[!valid] <- NA_real_
  if (!any(valid) || all(valid) || max_gap_frames == 0L) {
    return(list(values = values, status = factor(status, levels = FRAME_STATUS_LEVELS)))
  }

  observed_index <- which(valid)
  first_observed <- observed_index[1]
  last_observed <- observed_index[length(observed_index)]

  runs <- rle(valid)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  for (i in seq_along(runs$values)) {
    if (runs$values[i]) next
    gap_start <- starts[i]
    gap_end <- ends[i]
    # Interior gaps only: both neighbours must be observed.
    if (gap_start <= first_observed || gap_end >= last_observed) next
    if (runs$lengths[i] > max_gap_frames) next
    left <- gap_start - 1L
    right <- gap_end + 1L
    span <- right - left
    values[gap_start:gap_end] <- x[left] +
      (x[right] - x[left]) * (seq.int(gap_start, gap_end) - left) / span
    status[gap_start:gap_end] <- "interpolated"
  }

  list(values = values, status = factor(status, levels = FRAME_STATUS_LEVELS))
}

#' Apply bounded interpolation to tracked landmarks.
#'
#' x and y share one validity decision so a landmark is never half-observed.
#' Each landmark gains a `status` column recording observed / interpolated /
#' invalid, which `landmark_validity()` consumes.
#'
#' @param tracking a TrackingData object
#' @param landmarks landmarks to process; defaults to all
#' @param max_gap_s longest interior gap that may be interpolated, in seconds
#' @param likelihood_cutoff optional confidence below which a frame is rejected
#' @return the TrackingData object with filled coordinates and status columns
interpolate_tracking <- function(tracking, landmarks = names(tracking$data),
                                 max_gap_s = 0, likelihood_cutoff = NULL) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  fps <- get_tracking_fps(tracking)
  validate_scalar_number(max_gap_s, "max_gap_s", positive = TRUE, allow_zero = TRUE)
  if (!is.null(likelihood_cutoff)) {
    validate_scalar_number(likelihood_cutoff, "likelihood_cutoff")
  }
  max_gap_frames <- seconds_to_frames(max_gap_s, fps, round_fn = floor)

  for (point in landmarks) {
    values <- tracking$data[[point]]
    valid <- is.finite(values$x) & is.finite(values$y)
    if (!is.null(likelihood_cutoff) && !is.null(values$likelihood)) {
      confidence <- values$likelihood
      confidence[is.na(confidence)] <- -Inf
      valid <- valid & confidence >= likelihood_cutoff
    }
    filled_x <- bounded_linear_interpolation(values$x, valid, max_gap_frames)
    filled_y <- bounded_linear_interpolation(values$y, valid, max_gap_frames)
    tracking$data[[point]]$x <- filled_x$values
    tracking$data[[point]]$y <- filled_y$values
    # x and y are driven by the same validity mask, so their status agrees.
    tracking$data[[point]]$status <- filled_x$status
  }
  tracking$interpolation <- list(
    max_gap_s = max_gap_s,
    max_gap_frames = max_gap_frames,
    likelihood_cutoff = likelihood_cutoff,
    landmarks = landmarks
  )
  tracking
}

#' Per-landmark interpolation and missingness summary.
#'
#' @param tracking a TrackingData object processed by interpolate_tracking()
#' @return a data frame with one row per landmark
interpolation_report <- function(tracking, landmarks = names(tracking$data)) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  fps <- get_tracking_fps(tracking)
  do.call(rbind, lapply(landmarks, function(point) {
    status <- tracking$data[[point]]$status
    if (is.null(status)) {
      status <- factor(
        ifelse(
          is.finite(tracking$data[[point]]$x) & is.finite(tracking$data[[point]]$y),
          "observed", "invalid"
        ),
        levels = FRAME_STATUS_LEVELS
      )
    }
    invalid <- status == "invalid"
    gaps <- rle(as.logical(invalid))
    longest_gap <- if (any(gaps$values)) max(gaps$lengths[gaps$values]) else 0L
    n <- length(status)
    data.frame(
      landmark = point,
      frames = n,
      observed_frames = sum(status == "observed"),
      interpolated_frames = sum(status == "interpolated"),
      invalid_frames = sum(invalid),
      observed_fraction = sum(status == "observed") / n,
      interpolated_fraction = sum(status == "interpolated") / n,
      invalid_fraction = sum(invalid) / n,
      longest_invalid_gap_frames = as.integer(longest_gap),
      longest_invalid_gap_s = longest_gap / fps,
      stringsAsFactors = FALSE
    )
  }))
}
