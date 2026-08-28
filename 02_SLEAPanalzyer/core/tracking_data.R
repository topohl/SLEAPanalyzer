# Accessors and validation around the existing list-based TrackingData object.

validate_tracking_data <- function(tracking, required_landmarks = NULL) {
  if (!is.list(tracking) || !identical(tracking$object.type, "TrackingData")) {
    stop("Object is not a canonical TrackingData list")
  }
  if (!is.list(tracking$data) || length(tracking$data) == 0 ||
      is.null(names(tracking$data)) || any(!nzchar(names(tracking$data)))) {
    stop("TrackingData$data must be a named, non-empty list")
  }
  if (!is.numeric(tracking$frames) || length(tracking$frames) == 0 ||
      any(!is.finite(tracking$frames)) || anyDuplicated(tracking$frames)) {
    stop("TrackingData$frames must contain unique finite numeric frame identifiers")
  }
  validate_fps(tracking$fps)
  n_frames <- length(tracking$frames)

  for (point in names(tracking$data)) {
    coordinates <- tracking$data[[point]]
    if (!is.data.frame(coordinates) || !all(c("x", "y") %in% names(coordinates))) {
      stop("Tracked point '", point, "' must be a data frame with x and y columns")
    }
    if (!is.numeric(coordinates$x) || !is.numeric(coordinates$y) ||
        length(coordinates$x) != n_frames || length(coordinates$y) != n_frames) {
      stop("Tracked point '", point, "' is not aligned with TrackingData$frames")
    }
    if ("frame" %in% names(coordinates) &&
        !identical(as.numeric(coordinates$frame), as.numeric(tracking$frames))) {
      stop("Tracked point '", point, "' has frame identifiers that are not aligned")
    }
  }

  if (!is.null(tracking$seconds) && length(tracking$seconds) != n_frames) {
    stop("TrackingData$seconds is not aligned with TrackingData$frames")
  }
  if (!is.null(tracking$distance.units)) validate_coordinate_unit(tracking$distance.units)

  if (!is.null(required_landmarks)) {
    if (!is.character(required_landmarks) || anyNA(required_landmarks)) {
      stop("required_landmarks must be a character vector")
    }
    missing <- setdiff(required_landmarks, names(tracking$data))
    if (length(missing) > 0) {
      stop("TrackingData is missing landmark(s): ", paste(missing, collapse = ", "))
    }
  }
  invisible(TRUE)
}

get_tracking_frames <- function(tracking) {
  validate_tracking_data(tracking)
  tracking$frames
}

get_tracking_fps <- function(tracking) {
  validate_tracking_data(tracking)
  tracking$fps
}

get_tracking_duration <- function(tracking) {
  frames_to_seconds(length(get_tracking_frames(tracking)), get_tracking_fps(tracking))
}

get_tracking_unit <- function(tracking) {
  validate_tracking_data(tracking)
  if (is.null(tracking$distance.units)) {
    stop("TrackingData has no distance.units value")
  }
  tracking$distance.units
}

get_point_coordinates <- function(tracking, point, include_frame = FALSE,
                                  include_likelihood = FALSE) {
  validate_tracking_data(tracking, required_landmarks = point)
  columns <- c("x", "y")
  if (include_frame) columns <- c("frame", columns)
  if (include_likelihood && "likelihood" %in% names(tracking$data[[point]])) {
    columns <- c(columns, "likelihood")
  }
  tracking$data[[point]][, columns, drop = FALSE]
}

set_point_coordinates <- function(tracking, point, coordinates) {
  validate_tracking_data(tracking, required_landmarks = point)
  if (!is.data.frame(coordinates) || !all(c("x", "y") %in% names(coordinates))) {
    stop("coordinates must be a data frame with x and y columns")
  }
  n_frames <- length(get_tracking_frames(tracking))
  if (nrow(coordinates) != n_frames || !is.numeric(coordinates$x) || !is.numeric(coordinates$y)) {
    stop("coordinates must have one numeric x/y row per tracking frame")
  }
  tracking$data[[point]]$x <- coordinates$x
  tracking$data[[point]]$y <- coordinates$y
  if (!is.null(tracking$median.data) && point %in% rownames(tracking$median.data)) {
    tracking$median.data[point, "x"] <- stats::median(coordinates$x, na.rm = TRUE)
    tracking$median.data[point, "y"] <- stats::median(coordinates$y, na.rm = TRUE)
  }
  tracking
}

#' Per-frame validity mask for a set of landmarks.
#'
#' A frame is valid only if every requested landmark has finite coordinates on
#' that frame, and, when a likelihood cutoff is supplied, only if every
#' requested landmark meets it. Frames whose coordinates were interpolated are
#' reported as invalid when `interpolated_is_valid` is FALSE, which is the
#' conservative default for measurements that must not rest on fabricated data.
#'
#' @param tracking a TrackingData object
#' @param landmarks the landmarks a measurement depends on
#' @param likelihood_cutoff optional minimum tracking confidence
#' @param interpolated_is_valid whether interpolated frames count as observed
#' @return a logical vector, one entry per frame
landmark_validity <- function(tracking, landmarks,
                              likelihood_cutoff = NULL,
                              interpolated_is_valid = FALSE) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  if (!is.null(likelihood_cutoff)) {
    validate_scalar_number(likelihood_cutoff, "likelihood_cutoff")
  }
  n <- length(get_tracking_frames(tracking))
  valid <- rep(TRUE, n)
  for (point in landmarks) {
    values <- tracking$data[[point]]
    valid <- valid & is.finite(values$x) & is.finite(values$y)
    if (!is.null(likelihood_cutoff) && !is.null(values$likelihood)) {
      confidence <- values$likelihood
      confidence[is.na(confidence)] <- -Inf
      valid <- valid & confidence >= likelihood_cutoff
    }
    if (!isTRUE(interpolated_is_valid) && !is.null(values$status)) {
      valid <- valid & values$status != "interpolated"
    }
  }
  valid
}

has_landmarks <- function(tracking, landmarks) {
  validate_tracking_data(tracking)
  if (!is.character(landmarks) || anyNA(landmarks)) {
    stop("landmarks must be a character vector")
  }
  all(landmarks %in% names(tracking$data))
}
