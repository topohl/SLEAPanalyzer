# Repair of implausible landmarks.
#
# A pose model can put one landmark on the wrong thing for seconds at a time:
# the nose on an object outside the arena, an ear on the animal's reflection.
# The rest of the skeleton is then still right, and the misplaced point is
# inconsistent with it. This module finds such points by the skeleton itself,
# rebuilds the landmarks a measurement needs from rigid neighbours, and records
# what it did; bounded interpolation (interpolation.R) then bridges what is
# left. A rebuilt coordinate carries the status "reconstructed", so it is never
# counted as observed.

#' Landmarks of the body (not the tail, not the arena).
default_body_landmarks <- function() {
  c("nose", "leftEar", "rightEar", "spine1", "leftArm", "rightArm", "bodycentre",
    "leftFlank", "rightFlank", "spine2", "tailbase")
}

landmark_matrix <- function(tracking, landmarks) {
  lapply(stats::setNames(landmarks, landmarks), function(point) {
    xy <- get_point_coordinates(tracking, point)
    cbind(x = xy$x, y = xy$y)
  })
}

#' Landmarks inconsistent with the rest of the skeleton.
#'
#' For every pair of landmarks the typical distance is its median over the
#' recording. A pair is broken in a frame when its distance exceeds that median
#' by more than `max(margin, margin_fraction * median)`. A landmark is flagged
#' when it is broken against at least `min_share` of the other landmarks present
#' in that frame: a misplaced point disagrees with nearly every other point, a
#' correct one only with the misplaced ones. The medians are robust as long as
#' a landmark is misplaced in fewer than half of the frames.
#'
#' @param margin absolute tolerance, in coordinate units (cm after calibration)
#' @return logical matrix, frames x landmarks; attribute "medians" holds the
#'   typical pair distances
skeleton_outliers <- function(tracking, landmarks = default_body_landmarks(),
                              margin = 3, margin_fraction = 0.6, min_share = 0.5) {
  validate_tracking_data(tracking, required_landmarks = landmarks)
  validate_scalar_number(margin, "margin", positive = TRUE)
  validate_scalar_number(margin_fraction, "margin_fraction", positive = TRUE, allow_zero = TRUE)
  validate_scalar_number(min_share, "min_share", positive = TRUE)
  if (length(landmarks) < 3) stop("skeleton_outliers() needs at least three landmarks")
  L <- landmark_matrix(tracking, landmarks)
  pairs <- utils::combn(landmarks, 2)
  D <- vapply(seq_len(ncol(pairs)), function(k) {
    euclidean_distance(L[[pairs[1, k]]][, "x"], L[[pairs[1, k]]][, "y"], L[[pairs[2, k]]][, "x"], L[[pairs[2, k]]][, "y"])
  }, numeric(nrow(L[[1]])))
  D <- matrix(D, ncol = ncol(pairs))
  typical <- apply(D, 2, stats::median, na.rm = TRUE)
  limit <- typical + pmax(margin, margin_fraction * typical)
  known <- !is.na(D)
  broken <- known & sweep(D, 2, limit, ">")
  broken[is.na(broken)] <- FALSE
  flags <- vapply(landmarks, function(point) {
    k <- which(pairs[1, ] == point | pairs[2, ] == point)
    present <- rowSums(known[, k, drop = FALSE])
    present > 0 & rowSums(broken[, k, drop = FALSE]) >= ceiling(min_share * present)
  }, logical(nrow(D)))
  flags <- matrix(flags, ncol = length(landmarks), dimnames = list(NULL, landmarks))
  attr(flags, "medians") <- stats::setNames(typical, paste(pairs[1, ], pairs[2, ], sep = "-"))
  flags
}

#' Position of a landmark in the frame of two others.
#'
#' The target is expressed as `a + u * (b - a) + v * perp(b - a)`, with u and v
#' the medians over the frames where all three are trusted. `error_median` and
#' `error_p95` are the distances between prediction and observation over those
#' frames (coordinate units).
fit_anchor_model <- function(L, target, a, b, trusted) {
  ok <- trusted & stats::complete.cases(L[[target]], L[[a]], L[[b]])
  if (sum(ok) < 30) return(NULL)
  d <- L[[b]][ok, , drop = FALSE] - L[[a]][ok, , drop = FALSE]
  len2 <- rowSums(d^2); keep <- len2 > 0
  if (sum(keep) < 30) return(NULL)
  r <- L[[target]][ok, , drop = FALSE] - L[[a]][ok, , drop = FALSE]
  u <- (r[, 1] * d[, 1] + r[, 2] * d[, 2]) / len2
  v <- (-r[, 1] * d[, 2] + r[, 2] * d[, 1]) / len2
  u0 <- stats::median(u[keep]); v0 <- stats::median(v[keep])
  pred <- L[[a]][ok, , drop = FALSE] + u0 * d + v0 * cbind(-d[, 2], d[, 1])
  err <- sqrt(rowSums((pred - L[[target]][ok, , drop = FALSE])^2))
  data.frame(target = target, anchor_a = a, anchor_b = b, u = u0, v = v0, frames = sum(ok),
             error_median = stats::median(err), error_p95 = unname(stats::quantile(err, 0.95)),
             stringsAsFactors = FALSE)
}

predict_from_anchor <- function(L, model, rows) {
  a <- L[[model$anchor_a]][rows, , drop = FALSE]; b <- L[[model$anchor_b]][rows, , drop = FALSE]
  d <- b - a
  a + model$u * d + model$v * cbind(-d[, 2], d[, 1])
}

#' Find implausible landmarks, rebuild the ones a measurement needs, report.
#'
#' 1. skeleton_outliers() flags misplaced landmarks; their coordinates are
#'    removed. A frame with fewer than `min_landmarks` unflagged landmarks has
#'    a collapsed skeleton (the animal hidden or confused with something else):
#'    every body landmark of that frame is removed.
#' 2. Each `target` missing in a frame (removed above, or never tracked) is
#'    rebuilt from the best anchor pair present in that frame. Anchor models
#'    are fitted on the recording's own trusted frames, ranked by error, and
#'    only pairs whose 95th-percentile error is at most `max_anchor_error` are
#'    used. Rebuilt coordinates get status "reconstructed".
#' Run bounded interpolation afterwards for what remains missing.
#'
#' @return the TrackingData object with repaired coordinates; `$repair` holds
#'   the per-landmark summary, the episodes and the anchor models
repair_landmarks <- function(tracking, landmarks = default_body_landmarks(),
                             targets = c("nose", "bodycentre"), margin = 3,
                             margin_fraction = 0.6, min_share = 0.5,
                             min_landmarks = 4, max_anchor_error = 2) {
  validate_tracking_data(tracking, required_landmarks = unique(c(landmarks, targets)))
  if (!all(targets %in% landmarks)) stop("every repair target must be one of the checked landmarks")
  validate_positive_integer(min_landmarks, "min_landmarks")
  validate_scalar_number(max_anchor_error, "max_anchor_error", positive = TRUE)
  fps <- get_tracking_fps(tracking)
  flags <- skeleton_outliers(tracking, landmarks, margin, margin_fraction, min_share)
  L <- landmark_matrix(tracking, landmarks)
  present <- vapply(L, function(m) is.finite(m[, 1]) & is.finite(m[, 2]), logical(nrow(flags)))
  present <- matrix(present, ncol = length(landmarks), dimnames = dimnames(flags))
  collapse <- rowSums(present & !flags) < min_landmarks & rowSums(present) > 0
  removed <- (flags & present) | (collapse & present)
  trusted <- present & !removed
  status <- lapply(stats::setNames(landmarks, landmarks), function(point) {
    s <- tracking$data[[point]]$status
    if (is.null(s)) s <- ifelse(present[, point], "observed", "invalid")
    as.character(s)
  })
  for (point in landmarks) {
    gone <- removed[, point]
    if (any(gone)) {
      tracking$data[[point]]$x[gone] <- NA_real_; tracking$data[[point]]$y[gone] <- NA_real_
      status[[point]][gone] <- "invalid"
      L[[point]][gone, ] <- NA_real_
    }
  }
  models <- list(); rebuilt <- matrix(FALSE, nrow(flags), length(landmarks), dimnames = dimnames(flags))
  for (target in targets) {
    others <- setdiff(landmarks, target)
    cand <- utils::combn(others, 2)
    fits <- do.call(rbind, lapply(seq_len(ncol(cand)), function(k)
      fit_anchor_model(L, target, cand[1, k], cand[2, k], trusted[, target] & trusted[, cand[1, k]] & trusted[, cand[2, k]])))
    if (is.null(fits)) next
    fits <- fits[order(fits$error_p95), , drop = FALSE]
    usable <- fits[fits$error_p95 <= max_anchor_error, , drop = FALSE]
    models[[target]] <- fits
    need <- !is.finite(L[[target]][, 1]) & !collapse & rowSums(trusted) > 0
    for (k in seq_len(nrow(usable))) {
      if (!any(need)) break
      m <- usable[k, ]
      rows <- which(need & trusted[, m$anchor_a] & trusted[, m$anchor_b])
      if (!length(rows)) next
      xy <- predict_from_anchor(L, m, rows)
      tracking$data[[target]]$x[rows] <- xy[, 1]; tracking$data[[target]]$y[rows] <- xy[, 2]
      L[[target]][rows, ] <- xy
      status[[target]][rows] <- "reconstructed"
      rebuilt[rows, target] <- TRUE
      need[rows] <- FALSE
    }
  }
  for (point in landmarks) tracking$data[[point]]$status <- factor(status[[point]], levels = FRAME_STATUS_LEVELS)
  # flagged_s: the landmark alone was inconsistent; collapse_s: the whole skeleton was (the same for every landmark).
  summary <- data.frame(
    landmark = landmarks,
    flagged_s = colSums(flags & present & !collapse) / fps,
    collapse_s = colSums(collapse & present) / fps,
    reconstructed_s = colSums(rebuilt) / fps,
    still_missing_s = vapply(landmarks, function(p) sum(!is.finite(tracking$data[[p]]$x)) / fps, numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  episodes <- do.call(rbind, c(
    lapply(landmarks, function(point) {
      e <- episode_rows(removed[, point] & !collapse, fps, 0, check = "inconsistent_landmark", target = point)
      if (!nrow(e)) return(NULL)
      e$reconstructed_share <- vapply(seq_len(nrow(e)), function(i) mean(rebuilt[e$start_frame[i]:e$end_frame[i], point]), numeric(1))
      e }),
    list({ e <- episode_rows(collapse, fps, 0, check = "skeleton_collapse", target = "body")
           if (nrow(e)) { e$reconstructed_share <- 0; e } else NULL })))
  tracking$repair <- list(summary = summary, episodes = episodes, models = models,
                          frames_flagged = rowSums(flags & present) > 0, collapse = collapse,
                          settings = list(margin = margin, margin_fraction = margin_fraction, min_share = min_share,
                                          min_landmarks = min_landmarks, max_anchor_error = max_anchor_error))
  tracking
}
