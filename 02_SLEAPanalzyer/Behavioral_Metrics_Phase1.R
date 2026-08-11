# Small, dependency-free helpers shared by the Phase 1 NOR/SocP analyses and tests.

as_event_vector <- function(x) {
  x <- as.logical(x)
  x[is.na(x)] <- FALSE
  x
}

suppress_short_event_bouts <- function(x, min_frames = 1L) {
  x <- as_event_vector(x)
  if (length(min_frames) != 1 || !is.numeric(min_frames) ||
      !is.finite(min_frames) || min_frames < 1 || min_frames != floor(min_frames)) {
    stop("min_frames must be one positive integer")
  }
  if (length(x) == 0 || min_frames == 1) return(x)

  runs <- rle(x)
  runs$values <- runs$values & runs$lengths >= min_frames
  inverse.rle(runs)
}

event_entry_count <- function(x, min_frames = 1L) {
  x <- suppress_short_event_bouts(x, min_frames)
  if (length(x) == 0) return(0L)
  as.integer(sum(x & !c(FALSE, head(x, -1L))))
}

event_latency_s <- function(x, fps, min_frames = 1L) {
  if (length(fps) != 1 || !is.numeric(fps) || !is.finite(fps) || fps <= 0) {
    stop("fps must be one positive finite number")
  }
  x <- suppress_short_event_bouts(x, min_frames)
  first <- which(x)[1]
  if (is.na(first)) return(NA_real_)
  (first - 1) / fps
}

event_interbout_intervals_s <- function(x, fps, min_frames = 1L) {
  if (length(fps) != 1 || !is.numeric(fps) || !is.finite(fps) || fps <= 0) {
    stop("fps must be one positive finite number")
  }
  x <- suppress_short_event_bouts(x, min_frames)
  if (!any(x)) return(numeric())

  starts <- which(x & !c(FALSE, head(x, -1L)))
  ends <- which(x & !c(tail(x, -1L), FALSE))
  if (length(starts) < 2) return(numeric())
  (starts[-1L] - ends[-length(ends)] - 1L) / fps
}

event_summary <- function(x, fps, min_frames = 1L) {
  x <- suppress_short_event_bouts(x, min_frames)
  lengths <- if (any(x)) rle(x)$lengths[rle(x)$values] else numeric()
  list(
    duration_s = sum(x) / fps,
    bouts = event_entry_count(x),
    latency_s = event_latency_s(x, fps),
    mean_bout_s = if (length(lengths) == 0) NA_real_ else mean(lengths) / fps,
    max_bout_s = if (length(lengths) == 0) NA_real_ else max(lengths) / fps,
    interbout_intervals_s = event_interbout_intervals_s(x, fps)
  )
}

angle_between_vectors <- function(ax, ay, bx, by) {
  denominator <- sqrt(ax^2 + ay^2) * sqrt(bx^2 + by^2)
  cosine <- (ax * bx + ay * by) / denominator
  cosine[!is.finite(cosine)] <- NA_real_
  cosine <- pmax(-1, pmin(1, cosine))
  acos(cosine) * 180 / pi
}

tracking_xy <- function(tracking, point) {
  if (!point %in% names(tracking$data)) stop("Missing tracked point: ", point)
  data.frame(
    x = as.numeric(tracking$data[[point]]$x),
    y = as.numeric(tracking$data[[point]]$y)
  )
}

tracking_distance <- function(tracking, first, second) {
  a <- tracking_xy(tracking, first)
  b <- tracking_xy(tracking, second)
  if (nrow(a) != nrow(b)) stop("Tracked points have inconsistent frame counts")
  sqrt((a$x - b$x)^2 + (a$y - b$y)^2)
}

tracking_target_angle <- function(tracking, target) {
  nose <- tracking_xy(tracking, "nose")
  body <- tracking_xy(tracking, "bodycentre")
  target_xy <- tracking_xy(tracking, target)
  angle_between_vectors(
    nose$x - body$x,
    nose$y - body$y,
    nose$x - target_xy$x,
    nose$y - target_xy$y
  )
}

inside_axis_aligned_box <- function(x, y, center_x, center_y, width, height) {
  abs(x - center_x) <= width / 2 & abs(y - center_y) <= height / 2
}

validate_location_metadata <- function(location) {
  if (length(location) != 1 || is.na(location) || !location %in% c("L", "R")) {
    return(NA_character_)
  }
  as.character(location)
}

compute_nor_metrics <- function(tracking, novel_location, fps,
                                contact_distance = 4,
                                body_exclusion_distance = 1,
                                object_box_width = 9,
                                object_box_height = 7,
                                contact_angle = c(70, 290),
                                proximity_range = c(4, 8),
                                proximity_angle = c(90, 270)) {
  required <- c("nose", "bodycentre", "objL", "objR")
  missing <- setdiff(required, names(tracking$data))
  if (length(missing) > 0) stop("NOR tracking is missing: ", paste(missing, collapse = ", "))

  nose <- tracking_xy(tracking, "nose")
  obj_left <- tracking_xy(tracking, "objL")
  obj_right <- tracking_xy(tracking, "objR")
  n <- nrow(nose)
  if (any(c(nrow(obj_left), nrow(obj_right)) != n)) {
    stop("NOR tracked points have inconsistent frame counts")
  }

  location <- validate_location_metadata(novel_location)
  distance_left <- tracking_distance(tracking, "objL", "nose")
  distance_right <- tracking_distance(tracking, "objR", "nose")
  body_left <- tracking_distance(tracking, "objL", "bodycentre")
  body_right <- tracking_distance(tracking, "objR", "bodycentre")
  angle_left <- tracking_target_angle(tracking, "objL")
  angle_right <- tracking_target_angle(tracking, "objR")
  oriented_left <- abs(angle_left) >= contact_angle[1] & abs(angle_left) <= contact_angle[2]
  oriented_right <- abs(angle_right) >= contact_angle[1] & abs(angle_right) <= contact_angle[2]

  if (is.na(location)) {
    contact_left <- contact_right <- rep(NA, n)
  } else if (location == "R") {
    contact_left <- inside_axis_aligned_box(
      nose$x, nose$y, obj_left$x, obj_left$y, object_box_width, object_box_height
    ) & body_left > body_exclusion_distance & oriented_left
    contact_right <- distance_right <= contact_distance &
      body_right > body_exclusion_distance & oriented_right
  } else {
    contact_left <- distance_left <= contact_distance &
      body_left > body_exclusion_distance & oriented_left
    contact_right <- inside_axis_aligned_box(
      nose$x, nose$y, obj_right$x, obj_right$y, object_box_width, object_box_height
    ) & body_right > body_exclusion_distance & oriented_right
  }

  proximity_left <- distance_left > proximity_range[1] & distance_left <= proximity_range[2]
  proximity_right <- distance_right > proximity_range[1] & distance_right <= proximity_range[2]
  proximity_angle_left <- proximity_left &
    abs(angle_left) >= proximity_angle[1] & abs(angle_left) <= proximity_angle[2]
  proximity_angle_right <- proximity_right &
    abs(angle_right) >= proximity_angle[1] & abs(angle_right) <= proximity_angle[2]

  left <- if (is.na(location)) NULL else event_summary(contact_left, fps)
  right <- if (is.na(location)) NULL else event_summary(contact_right, fps)

  # Preserve the repository's legacy assignment: metadata R maps the left side to
  # novel and metadata L maps the right side to novel. Its biological meaning must
  # be confirmed before changing it.
  novel_is_left <- !is.na(location) && location == "R"
  mapped <- function(left_value, right_value) {
    if (is.na(location)) return(c(novel = NA_real_, familiar = NA_real_))
    if (novel_is_left) c(novel = left_value, familiar = right_value) else
      c(novel = right_value, familiar = left_value)
  }
  contact_mapped <- mapped(
    if (is.null(left)) NA_real_ else left$duration_s,
    if (is.null(right)) NA_real_ else right$duration_s
  )
  proximity_mapped <- mapped(sum(proximity_left, na.rm = TRUE) / fps,
                             sum(proximity_right, na.rm = TRUE) / fps)
  latency_values <- if (is.null(left)) numeric() else c(left$latency_s, right$latency_s)
  latency_values <- latency_values[is.finite(latency_values)]
  first_contact_latency <- if (length(latency_values) == 0) NA_real_ else min(latency_values)

  summary <- data.frame(
    contactLeft = if (is.null(left)) NA_real_ else left$duration_s,
    contactRight = if (is.null(right)) NA_real_ else right$duration_s,
    contactNov = unname(contact_mapped["novel"]),
    contactFam = unname(contact_mapped["familiar"]),
    proxLeft = sum(proximity_left, na.rm = TRUE) / fps,
    proxRight = sum(proximity_right, na.rm = TRUE) / fps,
    proxNov = unname(proximity_mapped["novel"]),
    proxFam = unname(proximity_mapped["familiar"]),
    proxLeftAngle = sum(proximity_angle_left, na.rm = TRUE) / fps,
    proxRightAngle = sum(proximity_angle_right, na.rm = TRUE) / fps,
    latency = first_contact_latency,
    latencyLeft = if (is.null(left)) NA_real_ else left$latency_s,
    latencyRight = if (is.null(right)) NA_real_ else right$latency_s,
    frequencyL = if (is.null(left)) NA_integer_ else left$bouts,
    frequencyR = if (is.null(right)) NA_integer_ else right$bouts,
    totalTime = n / fps,
    novelLoc = location
  )

  list(
    summary = summary,
    events = data.frame(
      contact_left = contact_left,
      contact_right = contact_right,
      proximity_left = proximity_left,
      proximity_right = proximity_right
    ),
    angles = data.frame(left = angle_left, right = angle_right),
    distances = data.frame(left = distance_left, right = distance_right)
  )
}

compute_socp_metrics <- function(tracking, novel_location, fps,
                                 contact_distance = 6,
                                 body_exclusion_distance = 1,
                                 proximity_range = c(6, 10)) {
  required <- c("nose", "bodycentre", "socl", "socr")
  missing <- setdiff(required, names(tracking$data))
  if (length(missing) > 0) stop("SocP tracking is missing: ", paste(missing, collapse = ", "))

  location <- validate_location_metadata(novel_location)
  distance_left <- tracking_distance(tracking, "socl", "nose")
  distance_right <- tracking_distance(tracking, "socr", "nose")
  body_left <- tracking_distance(tracking, "socl", "bodycentre")
  body_right <- tracking_distance(tracking, "socr", "bodycentre")
  contact_left <- distance_left <= contact_distance & body_left > body_exclusion_distance
  contact_right <- distance_right <= contact_distance & body_right > body_exclusion_distance
  proximity_left <- distance_left > proximity_range[1] & distance_left <= proximity_range[2]
  proximity_right <- distance_right > proximity_range[1] & distance_right <= proximity_range[2]
  left <- event_summary(contact_left, fps)
  right <- event_summary(contact_right, fps)

  novel_is_left <- !is.na(location) && location == "R"
  mapped <- function(left_value, right_value) {
    if (is.na(location)) return(c(novel = NA_real_, familiar = NA_real_))
    if (novel_is_left) c(novel = left_value, familiar = right_value) else
      c(novel = right_value, familiar = left_value)
  }
  contact_mapped <- mapped(left$duration_s, right$duration_s)
  proximity_mapped <- mapped(sum(proximity_left, na.rm = TRUE) / fps,
                             sum(proximity_right, na.rm = TRUE) / fps)
  latency_mapped <- mapped(left$latency_s, right$latency_s)

  list(
    summary = data.frame(
      contactNovel = unname(contact_mapped["novel"]),
      contactFamiliar = unname(contact_mapped["familiar"]),
      proxNovel = unname(proximity_mapped["novel"]),
      proxFamiliar = unname(proximity_mapped["familiar"]),
      contactLeft = left$duration_s,
      contactRight = right$duration_s,
      latencyNovel = unname(latency_mapped["novel"]),
      latencyFamiliar = unname(latency_mapped["familiar"]),
      latencyLeft = left$latency_s,
      latencyRight = right$latency_s,
      frequencyLeft = left$bouts,
      frequencyRight = right$bouts,
      totalTime = length(contact_left) / fps,
      novelLoc = location
    ),
    events = data.frame(
      contact_left = contact_left,
      contact_right = contact_right,
      proximity_left = proximity_left,
      proximity_right = proximity_right
    )
  )
}
