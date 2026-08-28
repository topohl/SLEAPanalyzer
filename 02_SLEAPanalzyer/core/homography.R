# Canonical arena rectification.
#
# Scalar x/y scaling is only correct when the camera looks straight down at the
# arena centre. Any tilt makes the pixels-per-centimetre ratio vary across the
# image, so a fixed scale factor over- or under-estimates distance depending on
# where in the arena the animal is. A projective homography removes that,
# because it is exactly the transform between a plane and its image under a
# pinhole camera.
#
# The output frame is defined so a rectangular arena spans (0, 0) to
# (width, height) in metric units, giving every assay one origin and one
# orientation regardless of how the camera was mounted.

#' Solve for the 3x3 homography mapping four source points to four targets.
#'
#' Uses the direct linear transform. With exactly four correspondences the
#' solution is exact; with more it is a least-squares fit.
#'
#' @param source a four-row data frame or matrix with x and y columns
#' @param target the corresponding destination points
#' @return a 3x3 matrix with H[3, 3] fixed to 1
solve_homography <- function(source, target) {
  source <- validate_correspondence(source, "source")
  target <- validate_correspondence(target, "target")
  if (nrow(source) != nrow(target)) {
    stop("source and target must have the same number of points")
  }
  if (nrow(source) < 4) stop("a homography requires at least four point pairs")

  n <- nrow(source)
  A <- matrix(0, nrow = 2 * n, ncol = 8)
  b <- numeric(2 * n)
  for (i in seq_len(n)) {
    x <- source$x[i]; y <- source$y[i]
    u <- target$x[i]; v <- target$y[i]
    A[2 * i - 1, ] <- c(x, y, 1, 0, 0, 0, -u * x, -u * y)
    A[2 * i, ]     <- c(0, 0, 0, x, y, 1, -v * x, -v * y)
    b[2 * i - 1] <- u
    b[2 * i] <- v
  }
  solution <- tryCatch(
    qr.solve(A, b),
    error = function(e) {
      stop("Arena corners are degenerate; cannot solve a homography: ",
           conditionMessage(e))
    }
  )
  matrix(c(solution, 1), nrow = 3, byrow = TRUE)
}

validate_correspondence <- function(points, name) {
  if (!is.data.frame(points) && !is.matrix(points)) {
    stop(name, " must be a data frame or matrix")
  }
  if (!all(c("x", "y") %in% colnames(points))) {
    stop(name, " must contain x and y columns")
  }
  out <- data.frame(x = as.numeric(points[, "x"]), y = as.numeric(points[, "y"]))
  if (any(!is.finite(out$x)) || any(!is.finite(out$y))) {
    stop(name, " must contain only finite coordinates")
  }
  out
}

#' Apply a homography to coordinate vectors.
#'
#' Points on or behind the camera plane (a non-positive homogeneous scale)
#' cannot be mapped and are returned as NA rather than as a reflected point.
apply_homography <- function(homography, x, y) {
  if (!is.matrix(homography) || !identical(dim(homography), c(3L, 3L))) {
    stop("homography must be a 3x3 matrix")
  }
  validate_numeric_vector(x, "x")
  validate_numeric_vector(y, "y")
  validate_equal_lengths(x, y, names = c("x", "y"))
  denominator <- homography[3, 1] * x + homography[3, 2] * y + homography[3, 3]
  denominator[!is.finite(denominator) | denominator <= 0] <- NA_real_
  list(
    x = (homography[1, 1] * x + homography[1, 2] * y + homography[1, 3]) / denominator,
    y = (homography[2, 1] * x + homography[2, 2] * y + homography[2, 3]) / denominator
  )
}

#' Canonical target corners for a rectangular arena.
#'
#' Corner order is top-left, top-right, bottom-right, bottom-left, matching the
#' repository convention of c("tl", "tr", "br", "bl"). The output frame has its
#' origin at the bottom-left corner with y increasing upward.
rectangular_arena_target <- function(width, height) {
  validate_scalar_number(width, "width", positive = TRUE)
  validate_scalar_number(height, "height", positive = TRUE)
  data.frame(
    x = c(0, width, width, 0),
    y = c(height, height, 0, 0)
  )
}

#' Build a calibration from four observed arena corners.
#'
#' @param corners observed corner coordinates in image units, ordered
#'   top-left, top-right, bottom-right, bottom-left
#' @param width,height the true arena dimensions in metric units
#' @param units the metric unit label, for example "cm"
#' @return an arena_calibration object carrying the transform and its QC
arena_calibration <- function(corners, width, height, units = "cm") {
  source <- validate_correspondence(corners, "corners")
  if (nrow(source) != 4) stop("a rectangular arena requires exactly four corners")
  canonical_coordinate_unit(units)
  target <- rectangular_arena_target(width, height)
  homography <- solve_homography(source, target)

  # Reprojection error: map the observed corners forward and compare against
  # where they should have landed. With four exact correspondences this is
  # numerically zero, so it detects solver failure rather than model misfit.
  projected <- apply_homography(homography, source$x, source$y)
  residuals <- sqrt((projected$x - target$x)^2 + (projected$y - target$y)^2)

  geometry <- rectangular_arena_geometry(source)
  # Opposite-side ratios reveal perspective: an untilted camera sees a
  # rectangle, so both ratios are 1.
  perspective <- max(abs(log(geometry$opposite_side_ratios)))

  structure(
    list(
      homography = homography,
      source_corners = source,
      target_corners = target,
      width = width,
      height = height,
      units = units,
      reprojection_error = residuals,
      max_reprojection_error = max(residuals),
      rms_reprojection_error = sqrt(mean(residuals^2)),
      geometry = geometry,
      perspective_index = perspective,
      valid = isTRUE(geometry$valid) && all(is.finite(residuals))
    ),
    class = "arena_calibration"
  )
}

#' Map image coordinates into canonical arena coordinates.
rectify_coordinates <- function(calibration, x, y) {
  if (!inherits(calibration, "arena_calibration")) {
    stop("calibration must be an arena_calibration object")
  }
  apply_homography(calibration$homography, x, y)
}

#' Rectify every landmark of a TrackingData object.
#'
#' Refuses to run twice, because the second application would map already
#' metric coordinates as though they were pixels.
#'
#' @param tracking a TrackingData object in image coordinates
#' @param calibration an arena_calibration object
#' @param landmarks landmarks to transform; defaults to all
rectify_tracking <- function(tracking, calibration, landmarks = names(tracking$data)) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  if (!inherits(calibration, "arena_calibration")) {
    stop("calibration must be an arena_calibration object")
  }
  if (!is.null(tracking$arena_calibration) ||
      (!is.null(tracking$distance.units) &&
       !identical(canonical_coordinate_unit(tracking$distance.units), "px"))) {
    stop("TrackingData is already calibrated; refusing to rectify twice")
  }
  if (!isTRUE(calibration$valid)) {
    stop("Arena calibration is not valid; refusing to rectify")
  }

  for (point in landmarks) {
    mapped <- rectify_coordinates(calibration, tracking$data[[point]]$x, tracking$data[[point]]$y)
    tracking$data[[point]]$x <- mapped$x
    tracking$data[[point]]$y <- mapped$y
  }
  tracking$distance.units <- calibration$units
  tracking$arena_calibration <- calibration
  if (!is.null(tracking$median.data)) {
    for (point in intersect(landmarks, rownames(tracking$median.data))) {
      tracking$median.data[point, "x"] <- stats::median(tracking$data[[point]]$x, na.rm = TRUE)
      tracking$median.data[point, "y"] <- stats::median(tracking$data[[point]]$y, na.rm = TRUE)
    }
  }
  tracking
}

#' Fraction of frames falling outside the calibrated arena.
#'
#' @param tolerance slack in metric units, so an animal pressed against a wall
#'   is not flagged by sub-centimetre landmark jitter
arena_violation_report <- function(tracking, landmarks = names(tracking$data),
                                   tolerance = 1) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  calibration <- tracking$arena_calibration
  if (is.null(calibration)) stop("TrackingData has no arena calibration")
  validate_scalar_number(tolerance, "tolerance", positive = TRUE, allow_zero = TRUE)

  do.call(rbind, lapply(landmarks, function(point) {
    x <- tracking$data[[point]]$x
    y <- tracking$data[[point]]$y
    observed <- is.finite(x) & is.finite(y)
    outside <- observed & (
      x < -tolerance | x > calibration$width + tolerance |
        y < -tolerance | y > calibration$height + tolerance
    )
    data.frame(
      landmark = point,
      observed_frames = sum(observed),
      outside_arena_frames = sum(outside),
      outside_arena_fraction = if (any(observed)) sum(outside) / sum(observed) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}
