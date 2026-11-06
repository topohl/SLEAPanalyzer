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

#' Distance of each point outside a polygon to its nearest edge.
#'
#' @param polygon matrix with columns x, y (the vertices in order)
#' @return 0 inside or on the polygon, NA where the point is unknown
distance_outside_polygon <- function(x, y, polygon) {
  d <- rep(NA_real_, length(x))
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(d)
  px <- x[ok]; py <- y[ok]; k <- nrow(polygon)
  edge <- vapply(seq_len(k), function(i) {
    a <- polygon[i, ]; b <- polygon[i %% k + 1, ]; dx <- b[1] - a[1]; dy <- b[2] - a[2]
    t <- pmin(1, pmax(0, ((px - a[1]) * dx + (py - a[2]) * dy) / (dx^2 + dy^2)))
    sqrt((px - a[1] - t * dx)^2 + (py - a[2] - t * dy)^2)
  }, numeric(length(px)))
  inside <- points_in_polygon(px, py, data.frame(x = polygon[, 1], y = polygon[, 2]))
  d[ok] <- ifelse(inside, 0, apply(matrix(edge, ncol = k), 1, min))
  d
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
#'   typical pair distances, "ratio" each pair's distance over its limit
#'   (frames x pairs, broken above 1) and "pairs" the pairs
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
  attr(flags, "ratio") <- sweep(D, 2, limit, "/")
  attr(flags, "pairs") <- pairs
  flags
}

#' A mutually consistent group of landmarks, frame by frame.
#'
#' Starting from `active`, the landmark broken against the most others still
#' in the group (ties: the larger summed excess over the limits) is removed,
#' one at a time, until no broken pair is left. When a skeleton is split
#' between the animal and something else, this keeps the larger part.
#'
#' @param ratio frames x pairs, distance over limit (attribute of skeleton_outliers())
#' @param active logical matrix, frames x landmarks: the landmarks to start from
#' @return logical matrix like `active`: the group kept
consistent_group <- function(ratio, pairs, landmarks, active) {
  ends <- vapply(landmarks, function(point) (pairs[1, ] == point | pairs[2, ] == point) * 1, numeric(ncol(pairs)))
  ends <- matrix(ends, ncol = length(landmarks))
  broken <- !is.na(ratio) & ratio > 1
  excess <- ifelse(broken, ratio - 1, 0)
  for (step in seq_along(landmarks)) {
    open <- broken & ((active * 1) %*% t(ends)) == 2
    rows <- which(rowSums(open) > 0)
    if (!length(rows)) break
    count <- (open[rows, , drop = FALSE] * 1) %*% ends
    extra <- (open[rows, , drop = FALSE] * excess[rows, , drop = FALSE]) %*% ends
    score <- count + extra / (extra + 1)
    score[!active[rows, , drop = FALSE]] <- -1
    active[cbind(rows, max.col(score, ties.method = "first"))] <- FALSE
  }
  active
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

predict_from_anchor <- function(L, model, rows, uv = NULL) {
  a <- L[[model$anchor_a]][rows, , drop = FALSE]; b <- L[[model$anchor_b]][rows, , drop = FALSE]
  d <- b - a
  if (is.null(uv)) uv <- cbind(rep(model$u, length(rows)), rep(model$v, length(rows)))
  a + uv[, 1] * d + uv[, 2] * cbind(-d[, 2], d[, 1])
}

#' The target's position in an anchor pair's frame of reference, carried into a gap.
#'
#' In frames where target and anchors were observed, and the anchors lie at
#' their usual distance apart (0.5 to 1.5 times its median), the target has a
#' position (u, v) relative to the anchor pair. For each frame of a gap, the
#' pose on either side is the median (u, v) of the `edge` nearest such frames;
#' the two sides are interpolated (or the one side that exists is taken) and
#' blended with the pair's median position: weight exp(-distance / decay) on
#' the carried pose, where distance is the time to the nearer such frame. Near
#' the edges of a gap the pose of the moment is kept; deep in a long gap the
#' typical pose takes over.
#' @param target0 the target's coordinates before any rebuild (missing where removed)
#' @param defined frames whose (u, v) may be measured
#' @param decay_frames decay of the weight, in frames
#' @param edge number of frames on each side whose median pose is carried
#' @return matrix of (u, v), one row per row of `rows`
pose_from_neighbours <- function(L, target0, model, rows, defined, decay_frames, edge = 5L) {
  uv0 <- cbind(rep(model$u, length(rows)), rep(model$v, length(rows)))
  if (decay_frames <= 0 || !length(rows)) return(uv0)
  d_all <- L[[model$anchor_b]] - L[[model$anchor_a]]; len <- sqrt(rowSums(d_all^2))
  typical <- stats::median(len[defined & is.finite(len) & len > 0])
  idx <- which(defined & is.finite(len) & len >= 0.5 * typical & len <= 1.5 * typical)
  if (!length(idx) || !is.finite(typical)) return(uv0)
  a <- L[[model$anchor_a]][idx, , drop = FALSE]; d <- d_all[idx, , drop = FALSE]
  r <- target0[idx, , drop = FALSE] - a; l2 <- rowSums(d^2)
  U <- cbind((r[, 1] * d[, 1] + r[, 2] * d[, 2]) / l2, (-r[, 1] * d[, 2] + r[, 2] * d[, 1]) / l2)
  n <- length(idx); k <- findInterval(rows, idx)                     # idx[k] < row < idx[k + 1]
  side <- function(lo, hi) {                                         # median pose of the defined frames lo..hi (indices into idx)
    out <- matrix(NA_real_, length(lo), 2); ok <- !is.na(lo) & !is.na(hi) & lo <= hi
    for (q in which(ok)) out[q, ] <- c(stats::median(U[lo[q]:hi[q], 1]), stats::median(U[lo[q]:hi[q], 2]))
    out
  }
  has_prev <- k >= 1; has_next <- k < n
  up <- side(ifelse(has_prev, pmax(1L, k - edge + 1L), NA), ifelse(has_prev, k, NA))
  un <- side(ifelse(has_next, k + 1L, NA), ifelse(has_next, pmin(n, k + edge), NA))
  dp <- ifelse(has_prev, rows - idx[pmax(k, 1L)], Inf); dn <- ifelse(has_next, idx[pmin(k + 1L, n)] - rows, Inf)
  both <- is.finite(dp) & is.finite(dn)
  frac <- ifelse(both, dp / (dp + dn), ifelse(is.finite(dp), 0, 1))
  up[!is.finite(dp), ] <- un[!is.finite(dp), ]; un[!is.finite(dn), ] <- up[!is.finite(dn), ]   # one side only: take it
  ui <- up * (1 - frac) + un * frac
  w <- exp(-pmin(dp, dn) / decay_frames); w[!is.finite(pmin(dp, dn))] <- 0
  out <- w * ui + (1 - w) * uv0
  bad <- !is.finite(out[, 1]) | !is.finite(out[, 2]); out[bad, ] <- uv0[bad, ]
  out
}

#' Find implausible landmarks, rebuild the ones a measurement needs, report.
#'
#' 0. With `corners`, a landmark more than `max_outside` outside the arena
#'    polygon is removed first: no part of the animal can be there.
#' 1. skeleton_outliers() flags misplaced landmarks; their coordinates are
#'    removed. A frame with fewer than `min_landmarks` unflagged landmarks
#'    keeps a consistent group that holds a majority of the landmarks present
#'    and at least `min_landmarks` (a split skeleton, consistent_group());
#'    otherwise the skeleton has collapsed (the
#'    animal hidden or confused with something else) and every body landmark
#'    of that frame is removed.
#' 2. Each `target` missing in a frame (removed above, or never tracked) is
#'    rebuilt from the best anchor pair present in that frame. Anchor models
#'    are fitted on the recording's own trusted frames, ranked by error, and
#'    only pairs whose 95th-percentile error is at most `max_anchor_error` are
#'    used. Rebuilt coordinates get status "reconstructed".
#'    With `head_check`, the nose is rebuilt only where both ears are present,
#'    and a rebuilt nose must sit where this animal's observed noses sit: its
#'    distance from the ear midpoint within their 1st to 99th percentile (as a
#'    share of the median head length) and its angle to the body axis (body
#'    centre to ear midpoint) at most their 99th percentile. A prediction that
#'    fails goes to the next anchor pair; a nose no pair can place stays missing
#'    (`refused_s`).
#'    With `pose_decay_s` > 0 and the head check, a rebuilt nose keeps the
#'    pose of the moment: its position relative to the anchor pair is carried
#'    in from the nearest frames where it was observed and fades to the pair's
#'    median position with time from those frames (pose_from_neighbours()). A
#'    nose whose pose of the moment fails the head check is placed at the
#'    pair's median position if that passes, before the next pair is tried.
#'    The body centre keeps the median position.
#' Run bounded interpolation afterwards for what remains missing.
#'
#' @return the TrackingData object with repaired coordinates; `$repair` holds
#'   the per-landmark summary, the episodes and the anchor models
repair_landmarks <- function(tracking, landmarks = default_body_landmarks(),
                             targets = c("nose", "bodycentre"), margin = 3,
                             margin_fraction = 0.6, min_share = 0.5,
                             min_landmarks = 4, max_anchor_error = 2,
                             corners = NULL, max_outside = NULL, head_check = TRUE, pose_decay_s = 1) {
  validate_tracking_data(tracking, required_landmarks = unique(c(landmarks, targets)))
  if (!all(targets %in% landmarks)) stop("every repair target must be one of the checked landmarks")
  validate_positive_integer(min_landmarks, "min_landmarks")
  validate_scalar_number(max_anchor_error, "max_anchor_error", positive = TRUE)
  validate_scalar_number(pose_decay_s, "pose_decay_s", positive = TRUE, allow_zero = TRUE)
  fps <- get_tracking_fps(tracking)
  n <- length(tracking$data[[landmarks[1]]]$x)
  far <- matrix(FALSE, n, length(landmarks), dimnames = list(NULL, landmarks))
  if (!is.null(corners)) {
    validate_scalar_number(max_outside, "max_outside", positive = TRUE)
    polygon <- do.call(rbind, lapply(corners, function(point) static_point(tracking, point)))
    for (point in landmarks) {
      d <- distance_outside_polygon(tracking$data[[point]]$x, tracking$data[[point]]$y, polygon)
      far[, point] <- !is.na(d) & d > max_outside
      tracking$data[[point]]$x[far[, point]] <- NA_real_; tracking$data[[point]]$y[far[, point]] <- NA_real_
    }
  }
  flags <- skeleton_outliers(tracking, landmarks, margin, margin_fraction, min_share)
  L <- landmark_matrix(tracking, landmarks)
  present <- vapply(L, function(m) is.finite(m[, 1]) & is.finite(m[, 2]), logical(nrow(flags)))
  present <- matrix(present, ncol = length(landmarks), dimnames = dimnames(flags))
  collapse <- rowSums(present & !flags) < min_landmarks & rowSums(present | far) > 0
  # When nearly half the landmarks are wrong, every landmark is broken against
  # half of the others. The frame is kept if one consistent group holds a
  # majority of the landmarks present, and at least min_landmarks (a split
  # skeleton); a few scattered points that happen to agree do not.
  split <- rep(FALSE, length(collapse))
  rows <- which(collapse)
  if (length(rows)) {
    kept <- consistent_group(attr(flags, "ratio")[rows, , drop = FALSE], attr(flags, "pairs"), landmarks,
                             present[rows, , drop = FALSE])
    clear <- rowSums(kept) >= min_landmarks & rowSums(kept) > rowSums(present[rows, , drop = FALSE]) / 2
    flags[rows[clear], ] <- present[rows[clear], , drop = FALSE] & !kept[clear, , drop = FALSE]
    split[rows[clear]] <- TRUE
    collapse[rows[clear]] <- FALSE
  }
  removed <- (flags & present) | (collapse & present) | far
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
  models <- list(); rebuilt <- matrix(FALSE, nrow(flags), length(landmarks), dimnames = dimnames(flags)); refused <- rebuilt
  head_ok <- isTRUE(head_check) && all(c("nose", "leftEar", "rightEar", "bodycentre") %in% landmarks)
  if (head_ok) {
    ears <- trusted[, "leftEar"] & trusted[, "rightEar"]
    em <- (L$leftEar + L$rightEar) / 2
    head_geometry <- function(xy) {   # distance from the ear midpoint and angle to the body axis
      h <- xy - em; ax <- em - L$bodycentre; len <- sqrt(rowSums(h^2))
      list(len = len, angle = acos(pmin(1, pmax(-1, rowSums(h * ax) / (len * sqrt(rowSums(ax^2)))))) * 180 / pi)
    }
    obs <- head_geometry(L$nose); ref <- trusted[, "nose"] & ears & is.finite(obs$len) & obs$len > 0
    head_ok <- sum(ref) >= 30
  }
  if (head_ok) {
    med <- stats::median(obs$len[ref]); len_lim <- stats::quantile(obs$len[ref] / med, c(.01, .99), names = FALSE)
    ang_max <- stats::quantile(obs$angle[ref & is.finite(obs$angle)], .99, names = FALSE)
  }
  for (target in targets) {
    target0 <- L[[target]]
    others <- setdiff(landmarks, target)
    cand <- utils::combn(others, 2)
    fits <- do.call(rbind, lapply(seq_len(ncol(cand)), function(k)
      fit_anchor_model(L, target, cand[1, k], cand[2, k], trusted[, target] & trusted[, cand[1, k]] & trusted[, cand[2, k]])))
    if (is.null(fits)) next
    fits <- fits[order(fits$error_p95), , drop = FALSE]
    usable <- fits[fits$error_p95 <= max_anchor_error, , drop = FALSE]
    models[[target]] <- fits
    need <- !is.finite(L[[target]][, 1]) & !collapse & rowSums(trusted) > 0
    check <- head_ok && target == "nose"
    if (check) { refused[need & !ears, target] <- TRUE; need <- need & ears }
    for (k in seq_len(nrow(usable))) {
      if (!any(need)) break
      m <- usable[k, ]
      rows <- which(need & trusted[, m$anchor_a] & trusted[, m$anchor_b])
      if (!length(rows)) next
      # the pose of the moment for the nose, which the head check guards; other targets keep the typical pose
      uv <- if (pose_decay_s > 0 && check) pose_from_neighbours(L, target0, m, rows, trusted[, target] & trusted[, m$anchor_a] & trusted[, m$anchor_b],
                                                                 pose_decay_s * fps) else NULL
      xy <- predict_from_anchor(L, m, rows, uv)
      finite <- is.finite(xy[, 1]) & is.finite(xy[, 2])
      if (!all(finite)) { rows <- rows[finite]; xy <- xy[finite, , drop = FALSE]; if (!is.null(uv)) uv <- uv[finite, , drop = FALSE]
        if (!length(rows)) next }
      if (check) {
        plausible <- function(xy) { full <- matrix(NA_real_, nrow(flags), 2); full[rows, ] <- xy; geo <- head_geometry(full)
          f <- geo$len[rows] / med >= len_lim[1] & geo$len[rows] / med <= len_lim[2] & (!is.finite(geo$angle[rows]) | geo$angle[rows] <= ang_max)
          f[is.na(f)] <- FALSE; f }
        fits <- plausible(xy)
        if (!is.null(uv) && any(!fits)) {
          # the pose of the moment can be implausible where the gap began (a head already half under a lid): the typical pose instead
          typical <- predict_from_anchor(L, m, rows)
          swap <- !fits & plausible(typical)
          xy[swap, ] <- typical[swap, , drop = FALSE]; fits <- fits | swap
        }
        refused[rows[!fits], target] <- TRUE
        rows <- rows[fits]; xy <- xy[fits, , drop = FALSE]
        if (!length(rows)) next
      }
      tracking$data[[target]]$x[rows] <- xy[, 1]; tracking$data[[target]]$y[rows] <- xy[, 2]
      L[[target]][rows, ] <- xy
      status[[target]][rows] <- "reconstructed"
      rebuilt[rows, target] <- TRUE
      need[rows] <- FALSE
    }
  }
  refused <- refused & !rebuilt
  for (point in landmarks) tracking$data[[point]]$status <- factor(status[[point]], levels = FRAME_STATUS_LEVELS)
  # outside_s: the landmark was outside the arena; flagged_s: it alone was inconsistent; collapse_s: the whole skeleton was
  # (frames in which the landmark was present; disjoint from outside_s).
  summary <- data.frame(
    landmark = landmarks,
    outside_s = colSums(far) / fps,
    flagged_s = colSums(flags & present & !collapse) / fps,
    collapse_s = colSums(collapse & present) / fps,
    reconstructed_s = colSums(rebuilt) / fps,
    still_missing_s = vapply(landmarks, function(p) sum(!is.finite(tracking$data[[p]]$x)) / fps, numeric(1)),
    refused_s = colSums(refused) / fps,
    stringsAsFactors = FALSE, row.names = NULL)
  landmark_episodes <- function(flag, check) lapply(landmarks, function(point) {
    e <- episode_rows(flag[, point], fps, 0, check = check, target = point)
    if (!nrow(e)) return(NULL)
    e$reconstructed_share <- vapply(seq_len(nrow(e)), function(i) mean(rebuilt[e$start_frame[i]:e$end_frame[i], point]), numeric(1))
    e })
  episodes <- do.call(rbind, c(
    landmark_episodes(far, "outside_arena"),
    landmark_episodes(removed & !far & !collapse, "inconsistent_landmark"),
    list({ e <- episode_rows(collapse, fps, 0, check = "skeleton_collapse", target = "body")
           if (nrow(e)) { e$reconstructed_share <- 0; e } else NULL },
         { e <- episode_rows(split, fps, 0, check = "skeleton_split", target = "body")
           if (nrow(e)) { e$reconstructed_share <- NA_real_; e } else NULL })))
  tracking$repair <- list(summary = summary, episodes = episodes, models = models,
                          frames_flagged = rowSums((flags & present) | far) > 0, collapse = collapse, split = split,
                          settings = list(margin = margin, margin_fraction = margin_fraction, min_share = min_share,
                                          min_landmarks = min_landmarks, max_anchor_error = max_anchor_error,
                                          max_outside = max_outside, head_check = head_check, pose_decay_s = pose_decay_s))
  tracking
}

#' Restore the length of a partly hidden head.
#'
#' An animal that pushes its head under an overhang (the lid of a wire cup)
#' keeps only the back of its head in view, and the pose model puts the nose
#' at the edge of what it can see: the tracked head is too short and the nose
#' sits short of where it is. Where the visible head (ear midpoint to nose) is
#' shorter than `below` times its median in the recording (frames whose nose
#' was observed), the nose is moved along the direction from the ear midpoint
#' to the nose until the head has that median length; frames without both
#' ears are left as they are. The nose moves by at most the median head
#' length. A moved nose has the status `reconstructed`. (A limit on the angle
#' between head and body axis was tried and removed: a mouse investigating a
#' cup turns its head, and the limit undid real recoveries.)
#'
#' Apply it after interpolation and after anything that measures the animal's
#' size (`apparent_body_size()` uses nose-ear separations), so that it changes
#' only measurements made from the nose.
#'
#' @param nose,ears the nose and the two ear landmarks
#' @param below a head shorter than this share of its median is restored
#' @return the TrackingData object with the nose moved; `$head_length` holds
#'   the median head length (`median`), the visible head length per frame
#'   (`length`), which frames were moved (`moved`) and the shift of the nose
#'   per frame (`shift`, 0 where unmoved)
restore_head_length <- function(tracking, nose = "nose", ears = c("leftEar", "rightEar"), below = 0.75) {
  validate_tracking_data(tracking, required_landmarks = c(nose, ears))
  n <- tracking$data[[nose]]; a <- tracking$data[[ears[1]]]; b <- tracking$data[[ears[2]]]
  ex <- (a$x + b$x) / 2; ey <- (a$y + b$y) / 2
  dx <- n$x - ex; dy <- n$y - ey; len <- sqrt(dx^2 + dy^2)
  status <- if (is.null(n$status)) rep("observed", length(len)) else as.character(n$status)
  reference <- is.finite(len) & len > 0 & status == "observed"
  none <- rep(FALSE, length(len))
  if (!any(reference)) {
    warning("restore_head_length: no frame with an observed nose and both ears; the nose is left as it is")
    tracking$head_length <- list(median = NA_real_, length = len, moved = none, shift = rep(0, length(len)))
    return(tracking)
  }
  m <- stats::median(len[reference])
  move <- is.finite(len) & len > 0 & len < below * m
  scale <- ifelse(move, m / len, 1)
  tracking$data[[nose]]$x <- ifelse(move, ex + dx * scale, n$x)
  tracking$data[[nose]]$y <- ifelse(move, ey + dy * scale, n$y)
  if (!is.null(n$status)) {
    s <- as.character(n$status); s[move] <- "reconstructed"
    tracking$data[[nose]]$status <- factor(s, levels = FRAME_STATUS_LEVELS)
  }
  tracking$head_length <- list(median = m, length = len, moved = move, shift = ifelse(move, m - len, 0))
  tracking
}
