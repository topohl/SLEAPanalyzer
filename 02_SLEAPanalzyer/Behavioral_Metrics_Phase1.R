# Compatibility and assay-definition layer for the shared behavioral core.

behavior_metrics_files <- unlist(lapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) character() else as.character(frame$ofile)
}), use.names = FALSE)
behavior_metrics_files <- behavior_metrics_files[
  basename(behavior_metrics_files) == "Behavioral_Metrics_Phase1.R"
]
behavior_metrics_root <- Sys.getenv("SLEAP_ANALYZER_REPO_ROOT")
behavior_metrics_dir <- if (length(behavior_metrics_files) > 0) {
  dirname(normalizePath(tail(behavior_metrics_files, 1), mustWork = TRUE))
} else if (nzchar(behavior_metrics_root)) {
  file.path(behavior_metrics_root, "02_SLEAPanalzyer")
} else {
  getwd()
}
source(file.path(behavior_metrics_dir, "core", "io.R"))
source_behavior_core(file.path(behavior_metrics_dir, "core"), envir = environment())
rm(behavior_metrics_files, behavior_metrics_root, behavior_metrics_dir)

# Phase 1 public names retained as delegating wrappers.
as_event_vector <- function(x) normalize_event_vector(x)

suppress_short_event_bouts <- function(x, min_frames = 1L) {
  filter_short_events(x, min_frames)
}

event_entry_count <- function(x, min_frames = 1L) {
  event_frequency(x, min_frames)
}

event_latency_s <- function(x, fps, min_frames = 1L) {
  event_latency(x, fps, min_frames)
}

event_interbout_intervals_s <- function(x, fps, min_frames = 1L) {
  event_interbout_intervals(x, fps, min_frames)
}

event_summary <- function(x, fps, min_frames = 1L) {
  summarize_event_metrics(x, fps, min_frames)
}

angle_between_vectors <- function(ax, ay, bx, by) {
  vector_angle_degrees(ax, ay, bx, by)
}

tracking_xy <- function(tracking, point) get_point_coordinates(tracking, point)

tracking_distance <- function(tracking, first, second) {
  tracking_point_distance(tracking, first, second)
}

tracking_target_angle <- function(tracking, target) {
  # Preserve the legacy NOR convention exactly: bodycentre -> nose compared
  # with target -> nose, rather than the more usual nose -> target vector.
  tracking_vector_angle(tracking, "bodycentre", "nose", target, "nose")
}

inside_axis_aligned_box <- function(x, y, center_x, center_y, width, height) {
  points_in_axis_aligned_box(x, y, center_x, center_y, width, height)
}

validate_location_metadata <- function(location) {
  if (length(location) != 1 || is.na(location) || !location %in% c("L", "R")) {
    return(NA_character_)
  }
  as.character(location)
}

validate_assay_timing <- function(tracking, fps) {
  validate_fps(fps)
  tracking_fps <- get_tracking_fps(tracking)
  if (!isTRUE(all.equal(as.numeric(fps), as.numeric(tracking_fps)))) {
    stop("fps does not match TrackingData$fps")
  }
  invisible(fps)
}

compute_nor_metrics <- function(tracking, novel_location, fps,
                                contact_distance = 4,
                                body_exclusion_distance = 1,
                                object_box_width = 9,
                                object_box_height = 7,
                                contact_angle = c(70, 290),
                                proximity_range = c(4, 8),
                                proximity_angle = c(90, 270),
                                threshold_unit = "cm") {
  required <- c("nose", "bodycentre", "objL", "objR")
  validate_tracking_data(tracking, required_landmarks = required)
  validate_assay_timing(tracking, fps)
  validate_threshold_unit(threshold_unit, get_tracking_unit(tracking))

  nose <- get_point_coordinates(tracking, "nose")
  obj_left <- get_point_coordinates(tracking, "objL")
  obj_right <- get_point_coordinates(tracking, "objR")
  n <- length(get_tracking_frames(tracking))
  location <- validate_location_metadata(novel_location)

  distance_left <- tracking_point_distance(tracking, "objL", "nose")
  distance_right <- tracking_point_distance(tracking, "objR", "nose")
  body_left <- tracking_point_distance(tracking, "objL", "bodycentre")
  body_right <- tracking_point_distance(tracking, "objR", "bodycentre")
  angle_left <- tracking_target_angle(tracking, "objL")
  angle_right <- tracking_target_angle(tracking, "objR")
  oriented_left <- abs(angle_left) >= contact_angle[1] & abs(angle_left) <= contact_angle[2]
  oriented_right <- abs(angle_right) >= contact_angle[1] & abs(angle_right) <= contact_angle[2]

  if (is.na(location)) {
    contact_left <- contact_right <- rep(NA, n)
  } else if (location == "R") {
    contact_left <- points_in_axis_aligned_box(
      nose$x, nose$y, obj_left$x, obj_left$y, object_box_width, object_box_height
    ) & body_left > body_exclusion_distance & oriented_left
    contact_right <- distance_right <= contact_distance &
      body_right > body_exclusion_distance & oriented_right
  } else {
    contact_left <- distance_left <= contact_distance &
      body_left > body_exclusion_distance & oriented_left
    contact_right <- points_in_axis_aligned_box(
      nose$x, nose$y, obj_right$x, obj_right$y, object_box_width, object_box_height
    ) & body_right > body_exclusion_distance & oriented_right
  }

  proximity_left <- distance_left > proximity_range[1] & distance_left <= proximity_range[2]
  proximity_right <- distance_right > proximity_range[1] & distance_right <= proximity_range[2]
  proximity_angle_left <- proximity_left &
    abs(angle_left) >= proximity_angle[1] & abs(angle_left) <= proximity_angle[2]
  proximity_angle_right <- proximity_right &
    abs(angle_right) >= proximity_angle[1] & abs(angle_right) <= proximity_angle[2]

  left <- if (is.na(location)) NULL else summarize_event_metrics(contact_left, fps)
  right <- if (is.na(location)) NULL else summarize_event_metrics(contact_right, fps)

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
  proximity_mapped <- mapped(
    frames_to_seconds(sum(proximity_left, na.rm = TRUE), fps),
    frames_to_seconds(sum(proximity_right, na.rm = TRUE), fps)
  )
  latency_values <- if (is.null(left)) numeric() else c(left$latency_s, right$latency_s)
  latency_values <- latency_values[is.finite(latency_values)]
  first_contact_latency <- if (length(latency_values) == 0) NA_real_ else min(latency_values)

  summary <- data.frame(
    contactLeft = if (is.null(left)) NA_real_ else left$duration_s,
    contactRight = if (is.null(right)) NA_real_ else right$duration_s,
    contactNov = unname(contact_mapped["novel"]),
    contactFam = unname(contact_mapped["familiar"]),
    proxLeft = frames_to_seconds(sum(proximity_left, na.rm = TRUE), fps),
    proxRight = frames_to_seconds(sum(proximity_right, na.rm = TRUE), fps),
    proxNov = unname(proximity_mapped["novel"]),
    proxFam = unname(proximity_mapped["familiar"]),
    proxLeftAngle = frames_to_seconds(sum(proximity_angle_left, na.rm = TRUE), fps),
    proxRightAngle = frames_to_seconds(sum(proximity_angle_right, na.rm = TRUE), fps),
    latency = first_contact_latency,
    latencyLeft = if (is.null(left)) NA_real_ else left$latency_s,
    latencyRight = if (is.null(right)) NA_real_ else right$latency_s,
    frequencyL = if (is.null(left)) NA_integer_ else left$bouts,
    frequencyR = if (is.null(right)) NA_integer_ else right$bouts,
    totalTime = get_tracking_duration(tracking),
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
