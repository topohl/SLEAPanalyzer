# Frame-level plausibility checks.
#
# Session QC (qc.R) decides whether a recording is usable. The checks here find
# the stretches of a usable recording that the scoring rules do not expect: the
# animal raised off the floor onto an object, the body in a place a floor-level
# animal cannot reach, implausible jumps, a distorted skeleton. Each check
# returns per-frame flags and an episode table. Like qc.R they never exclude
# anything by themselves; an assay applies a check only when its configuration
# says so, and records what it applied.

# ---------------------------------------------------------------------------
# Apparent body size
# ---------------------------------------------------------------------------

#' Landmark pairs whose separation changes little with posture.
default_size_pairs <- function() {
  list(
    c("leftEar", "rightEar"), c("nose", "leftEar"), c("nose", "rightEar"),
    c("leftFlank", "rightFlank"), c("spine1", "bodycentre"), c("bodycentre", "spine2")
  )
}

#' Apparent body size relative to the recording's own median.
#'
#' A top-down camera images an animal that is closer to the lens larger. Each
#' landmark pair's separation is divided by its median over the recording, and
#' a frame's value is the median of those ratios, so one misplaced landmark
#' cannot carry it. On the floor the value is about 1; a mouse on the lid of a
#' wire cup in the Exp9 social-preference arena measures about 1.15 to 1.25.
#'
#' @param pairs list of length-two landmark names; pairs naming a landmark the
#'   tracking data lack are skipped
#' @param min_pairs frames with fewer measurable pairs are NA
#' @return numeric vector, one value per frame
apparent_body_size <- function(tracking, pairs = default_size_pairs(), min_pairs = 2L) {
  validate_tracking_data(tracking)
  validate_positive_integer(min_pairs, "min_pairs")
  pairs <- Filter(function(pair) all(pair %in% names(tracking$data)), pairs)
  if (length(pairs) < min_pairs) {
    stop("apparent_body_size() needs at least ", min_pairs,
         " landmark pairs present in the tracking data")
  }
  n <- length(get_tracking_frames(tracking))
  ratios <- vapply(pairs, function(pair) {
    separation <- tracking_point_distance(tracking, pair[1], pair[2])
    typical <- stats::median(separation[is.finite(separation) & separation > 0])
    if (!is.finite(typical)) rep(NA_real_, n) else separation / typical
  }, numeric(n))
  ratios <- matrix(ratios, nrow = n)
  measurable <- rowSums(is.finite(ratios))
  size <- apply(ratios, 1, stats::median, na.rm = TRUE)
  size[measurable < min_pairs] <- NA_real_
  size
}

# ---------------------------------------------------------------------------
# Position relative to a static object
# ---------------------------------------------------------------------------

static_point <- function(tracking, point) {
  coordinates <- get_point_coordinates(tracking, point)
  c(x = stats::median(coordinates$x, na.rm = TRUE), y = stats::median(coordinates$y, na.rm = TRUE))
}

#' Body position relative to a static object.
#'
#' `distance` is the distance from the object point. `beyond` is the component
#' along the direction from the object to the nearest arena corner: positive
#' values lie between the object and the walls of that corner.
#'
#' @param object landmark name of the object (its median position is used)
#' @param corners landmark names of the arena corners
#' @param body landmark used as the body position
object_relative_position <- function(tracking, object, corners, body = "bodycentre") {
  validate_tracking_data(tracking, required_landmarks = c(object, corners, body))
  centre <- static_point(tracking, object)
  corner_xy <- do.call(rbind, lapply(corners, function(point) static_point(tracking, point)))
  nearest <- which.min((corner_xy[, "x"] - centre["x"])^2 + (corner_xy[, "y"] - centre["y"])^2)
  direction <- corner_xy[nearest, ] - centre
  direction <- direction / sqrt(sum(direction^2))
  position <- get_point_coordinates(tracking, body)
  dx <- position$x - centre["x"]
  dy <- position$y - centre["y"]
  data.frame(
    distance = sqrt(dx^2 + dy^2),
    beyond = dx * direction["x"] + dy * direction["y"],
    row.names = NULL
  )
}

# ---------------------------------------------------------------------------
# On-object detection
# ---------------------------------------------------------------------------

empty_object_episode_table <- function() {
  data.frame(
    object = character(), start_frame = integer(), end_frame = integer(),
    start_s = numeric(), duration_s = numeric(), candidate_frames = integer(),
    median_apparent_size = numeric(), min_distance = numeric(), max_beyond = numeric(),
    decision = character(), stringsAsFactors = FALSE
  )
}

#' Running median that ignores missing values.
rolling_median <- function(x, window_frames) {
  validate_positive_integer(window_frames, "window_frames")
  n <- length(x)
  half <- window_frames %/% 2L
  vapply(seq_len(n), function(i) {
    v <- x[max(1L, i - half):min(n, i + half)]
    v <- v[is.finite(v)]
    if (length(v) == 0) NA_real_ else stats::median(v)
  }, numeric(1))
}

#' Episodes in which the animal is on top of a static object.
#'
#' Candidate frames put the body where a floor-level animal cannot be: within
#' `footprint_radius` of the object point (over the object) or, when `beyond`
#' is given, more than `beyond` past the object toward the nearest arena corner
#' (between the object and the walls), always within `max_distance`. Candidate
#' runs separated by at most `merge_gap_s` form one episode. An episode is
#' accepted as on the object when it lasts at least `min_episode_s` and the
#' animal looks raised: the episode's median apparent body size is at least
#' `min_apparent_size`. An animal wedged beside an object on the floor, or a
#' tracking glitch that drops the body point onto the object, lands in the same
#' places at normal size, so those episodes are reported but not accepted.
#' An accepted episode is then extended backward and forward for as long as the
#' animal stays within `max_distance` and still looks raised (apparent size,
#' as a running median over `size_window_s`, at least `min_apparent_size`):
#' this covers the lid edge, where the body is not over the object.
#'
#' @param size optional apparent_body_size() series, to avoid recomputing it
#'   for several objects
#' @return list: `on_object`, a logical vector (every accepted episode,
#'   extended), and `episodes`, one row per candidate episode with its decision
#'   ("on_object", "normal_size" or "too_short"); `start_frame`/`end_frame`
#'   of accepted episodes include the extension
detect_on_object <- function(tracking, object, corners, footprint_radius,
                             beyond = NULL, max_distance = Inf,
                             min_apparent_size = 1.08, merge_gap_s = 0.5,
                             min_episode_s = 0.2, size_window_s = 0.5,
                             body = "bodycentre", size = NULL) {
  validate_scalar_number(footprint_radius, "footprint_radius", positive = TRUE)
  if (!is.null(beyond)) validate_scalar_number(beyond, "beyond", allow_zero = TRUE)
  validate_scalar_number(min_apparent_size, "min_apparent_size", positive = TRUE)
  validate_scalar_number(merge_gap_s, "merge_gap_s", positive = TRUE, allow_zero = TRUE)
  validate_scalar_number(min_episode_s, "min_episode_s", positive = TRUE, allow_zero = TRUE)
  fps <- get_tracking_fps(tracking)
  n <- length(get_tracking_frames(tracking))
  if (is.null(size)) size <- apparent_body_size(tracking)
  validate_equal_lengths(size, seq_len(n), names = c("size", "frames"))

  position <- object_relative_position(tracking, object, corners, body)
  candidate <- position$distance < footprint_radius
  if (!is.null(beyond)) candidate <- candidate | position$beyond > beyond
  candidate <- candidate & position$distance < max_distance
  candidate[is.na(candidate)] <- FALSE

  on_object <- rep(FALSE, n)
  merged <- merge_event_gaps(candidate, seconds_to_frames(merge_gap_s, fps, round_fn = floor))
  bouts <- event_bout_table(merged, fps)
  if (nrow(bouts) == 0) {
    return(list(on_object = on_object, episodes = empty_object_episode_table()))
  }
  min_frames <- seconds_to_frames(min_episode_s, fps)
  episodes <- do.call(rbind, lapply(seq_len(nrow(bouts)), function(i) {
    span <- bouts$start_frame[i]:bouts$end_frame[i]
    inside <- span[candidate[span]]
    typical_size <- stats::median(size[inside], na.rm = TRUE)
    decision <- if (length(span) < min_frames) "too_short" else
      if (is.finite(typical_size) && typical_size >= min_apparent_size) "on_object" else "normal_size"
    data.frame(
      object = object, start_frame = bouts$start_frame[i], end_frame = bouts$end_frame[i],
      start_s = bouts$latency_seconds[i], duration_s = bouts$duration_seconds[i],
      candidate_frames = length(inside), median_apparent_size = typical_size,
      min_distance = min(position$distance[inside], na.rm = TRUE),
      max_beyond = max(position$beyond[inside], na.rm = TRUE),
      decision = decision, stringsAsFactors = FALSE
    )
  }))
  accepted <- which(episodes$decision == "on_object")
  if (length(accepted) > 0) {
    raised <- rolling_median(size, max(1L, seconds_to_frames(size_window_s, fps, round_fn = round)))
    raised <- is.finite(raised) & raised >= min_apparent_size &
      is.finite(position$distance) & position$distance < max_distance
    for (i in accepted) {
      first <- episodes$start_frame[i]
      last <- episodes$end_frame[i]
      while (first > 1L && raised[first - 1L]) first <- first - 1L
      while (last < n && raised[last + 1L]) last <- last + 1L
      episodes$start_frame[i] <- first
      episodes$end_frame[i] <- last
      episodes$start_s[i] <- (first - 1L) / fps
      episodes$duration_s[i] <- (last - first + 1L) / fps
      on_object[first:last] <- TRUE
    }
  }
  list(on_object = on_object, episodes = episodes)
}

# ---------------------------------------------------------------------------
# Reviewed decisions
# ---------------------------------------------------------------------------

#' Read reviewed QC decisions.
#'
#' One row per reviewed stretch: `file` (the recording name without
#' extension), `target` (the object landmark), `start_s`, `end_s` and
#' `decision`: "on_object" forces the stretch to count as on the object,
#' "not_on_object" forces it not to. `reason` and `reviewer` are kept for the
#' record. An absent file means no decisions.
read_qc_decisions <- function(path) {
  columns <- c("file", "target", "start_s", "end_s", "decision", "reason", "reviewer")
  if (is.null(path) || !file.exists(path)) {
    return(stats::setNames(data.frame(matrix(nrow = 0, ncol = length(columns))), columns))
  }
  decisions <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(columns, names(decisions))
  if (length(missing) > 0) stop("QC decisions file lacks column(s): ", paste(missing, collapse = ", "))
  bad <- setdiff(unique(decisions$decision), c("on_object", "not_on_object"))
  if (length(bad) > 0) stop("Unknown QC decision(s): ", paste(bad, collapse = ", "))
  if (any(!is.finite(decisions$start_s) | !is.finite(decisions$end_s) | decisions$end_s < decisions$start_s)) {
    stop("QC decisions need finite start_s <= end_s")
  }
  decisions[, columns]
}

#' Apply reviewed decisions to an on-object mask.
#'
#' Frame i (1-based) lies at (i - 1) / fps seconds, the convention of every
#' latency in this repository.
apply_qc_decisions <- function(on_object, decisions, file, target, fps) {
  rows <- decisions[decisions$file == file & decisions$target == target, , drop = FALSE]
  if (nrow(rows) == 0) return(on_object)
  time_s <- (seq_along(on_object) - 1L) / fps
  for (i in seq_len(nrow(rows))) {
    span <- time_s >= rows$start_s[i] & time_s <= rows$end_s[i]
    on_object[span] <- rows$decision[i] == "on_object"
  }
  on_object
}

# ---------------------------------------------------------------------------
# Generic report-only checks
# ---------------------------------------------------------------------------

#' Episodes of implausible displacement of one landmark (jumps).
#'
#' @param max_speed maximum plausible speed in coordinate units per second
jump_episodes <- function(tracking, landmark = "bodycentre", max_speed, merge_gap_s = 0.2) {
  validate_scalar_number(max_speed, "max_speed", positive = TRUE)
  fps <- get_tracking_fps(tracking)
  xy <- get_point_coordinates(tracking, landmark)
  speed <- c(NA_real_, euclidean_distance(xy$x[-1], xy$y[-1], head(xy$x, -1), head(xy$y, -1))) * fps
  flag <- is.finite(speed) & speed > max_speed
  episode_rows(flag, fps, merge_gap_s, check = "jump", target = landmark,
               statistic = speed, statistic_name = "max_speed", summary = max)
}

#' Episodes in which the body axis is far longer or shorter than usual.
skeleton_episodes <- function(tracking, front = "nose", back = "tailbase",
                              tolerance = 0.5, merge_gap_s = 0.2, min_episode_s = 0.2) {
  fps <- get_tracking_fps(tracking)
  length_series <- tracking_point_distance(tracking, front, back)
  ratio <- length_series / stats::median(length_series[is.finite(length_series)])
  flag <- is.finite(ratio) & abs(ratio - 1) > tolerance
  rows <- episode_rows(flag, fps, merge_gap_s, check = "skeleton",
                       target = paste(front, back, sep = "-"),
                       statistic = ratio, statistic_name = "median_length_ratio",
                       summary = stats::median)
  rows[rows$duration_s >= min_episode_s, , drop = FALSE]
}

#' Episodes in which a landmark lies outside the arena polygon.
#'
#' @param tolerance distance outside the polygon tolerated, in coordinate units
outside_arena_episodes <- function(tracking, landmark = "bodycentre", corners,
                                   tolerance = 0, merge_gap_s = 0.2, min_episode_s = 0.2) {
  fps <- get_tracking_fps(tracking)
  polygon <- do.call(rbind, lapply(corners, function(point) static_point(tracking, point)))
  if (tolerance > 0) {
    centre <- colMeans(polygon)
    scale <- 1 + tolerance / mean(sqrt(rowSums(sweep(polygon, 2, centre)^2)))
    polygon <- sweep(sweep(polygon, 2, centre) * scale, 2, centre, "+")
  }
  xy <- get_point_coordinates(tracking, landmark)
  inside <- points_in_polygon(xy$x, xy$y, data.frame(x = polygon[, "x"], y = polygon[, "y"]))
  flag <- is.finite(xy$x) & is.finite(xy$y) & !inside
  rows <- episode_rows(flag, fps, merge_gap_s, check = "outside_arena", target = landmark)
  rows[rows$duration_s >= min_episode_s, , drop = FALSE]
}

#' One row per episode of a flag vector, in the common QC episode layout.
episode_rows <- function(flag, fps, merge_gap_s, check, target,
                         statistic = NULL, statistic_name = NA_character_, summary = NULL) {
  merged <- merge_event_gaps(flag, seconds_to_frames(merge_gap_s, fps, round_fn = floor))
  bouts <- event_bout_table(merged, fps)
  if (nrow(bouts) == 0) {
    return(data.frame(check = character(), target = character(), start_frame = integer(),
                      end_frame = integer(), start_s = numeric(), duration_s = numeric(),
                      statistic_name = character(), statistic = numeric(), stringsAsFactors = FALSE))
  }
  value <- if (is.null(statistic)) rep(NA_real_, nrow(bouts)) else vapply(seq_len(nrow(bouts)), function(i) {
    span <- bouts$start_frame[i]:bouts$end_frame[i]
    v <- statistic[span][flag[span] & is.finite(statistic[span])]
    if (length(v) == 0) NA_real_ else summary(v)
  }, numeric(1))
  data.frame(
    check = check, target = target, start_frame = bouts$start_frame, end_frame = bouts$end_frame,
    start_s = bouts$latency_seconds, duration_s = bouts$duration_seconds,
    statistic_name = statistic_name, statistic = value, stringsAsFactors = FALSE
  )
}
