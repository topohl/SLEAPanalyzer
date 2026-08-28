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

#' Score nose contact with one object.
#'
#' The same detector must be applied to both objects. Using a different
#' detection region for the novel and the familiar object biases the
#' discrimination index by construction, because the regions have different
#' areas.
#'
#' @param nose nose coordinates
#' @param object object coordinates
#' @param geometry "radial" (nose within `radius` of the object point) or
#'   "box" (nose inside a `width` x `height` axis-aligned footprint)
#' @return a logical vector, NA where coordinates are missing
object_contact_region <- function(nose, object, geometry = c("radial", "box"),
                                  radius = NULL, width = NULL, height = NULL) {
  geometry <- match.arg(geometry)
  if (geometry == "radial") {
    validate_scalar_number(radius, "contact_distance", positive = TRUE)
    euclidean_distance(nose$x, nose$y, object$x, object$y) <= radius
  } else {
    validate_scalar_number(width, "object_box_width", positive = TRUE)
    validate_scalar_number(height, "object_box_height", positive = TRUE)
    points_in_axis_aligned_box(nose$x, nose$y, object$x, object$y, width, height)
  }
}

#' Novel-object-recognition measurements.
#'
#' @param contact_geometry Detection region applied to *both* objects.
#'   "radial" uses `contact_distance`; "box" uses `object_box_width` /
#'   `object_box_height`. "legacy_asymmetric" reproduces the biased pre-v2
#'   behavior and is provided only for reanalysis of existing outputs.
#' @param valid optional per-frame validity mask; defaults to frames where
#'   every required landmark has finite coordinates.
#' @param min_bout_s minimum contact bout duration, in seconds.
#' @param max_gap_s interior gap bridged before bouts are filtered, in seconds.
compute_nor_metrics <- function(tracking, novel_location, fps,
                                contact_distance = 4,
                                body_exclusion_distance = 1,
                                object_box_width = 9,
                                object_box_height = 7,
                                contact_angle = c(70, 290),
                                proximity_range = c(4, 8),
                                proximity_angle = c(90, 270),
                                threshold_unit = "cm",
                                contact_geometry = c("radial", "box", "legacy_asymmetric"),
                                valid = NULL,
                                min_bout_s = 0,
                                max_gap_s = 0) {
  required <- c("nose", "bodycentre", "objL", "objR")
  validate_tracking_data(tracking, required_landmarks = required)
  validate_assay_timing(tracking, fps)
  validate_threshold_unit(threshold_unit, get_tracking_unit(tracking))
  contact_geometry <- match.arg(contact_geometry)

  nose <- get_point_coordinates(tracking, "nose")
  obj_left <- get_point_coordinates(tracking, "objL")
  obj_right <- get_point_coordinates(tracking, "objR")
  n <- length(get_tracking_frames(tracking))
  location <- validate_location_metadata(novel_location)
  if (is.null(valid)) valid <- landmark_validity(tracking, required)
  valid <- normalize_validity_mask(valid, n)

  distance_left <- tracking_point_distance(tracking, "objL", "nose")
  distance_right <- tracking_point_distance(tracking, "objR", "nose")
  body_left <- tracking_point_distance(tracking, "objL", "bodycentre")
  body_right <- tracking_point_distance(tracking, "objR", "bodycentre")
  angle_left <- tracking_target_angle(tracking, "objL")
  angle_right <- tracking_target_angle(tracking, "objR")
  oriented_left <- abs(angle_left) >= contact_angle[1] & abs(angle_left) <= contact_angle[2]
  oriented_right <- abs(angle_right) >= contact_angle[1] & abs(angle_right) <= contact_angle[2]

  if (contact_geometry == "legacy_asymmetric") {
    warning(
      "contact_geometry = 'legacy_asymmetric' scores the novel object with a ",
      object_box_width, "x", object_box_height, " box and the familiar object ",
      "with a ", contact_distance, " radius. The regions have different areas, ",
      "so the discrimination index is biased. Use it only to reproduce ",
      "pre-v2 outputs."
    )
    if (is.na(location)) {
      region_left <- region_right <- rep(NA, n)
    } else {
      box_side <- if (location == "R") "left" else "right"
      region_left <- if (box_side == "left") {
        object_contact_region(nose, obj_left, "box",
                              width = object_box_width, height = object_box_height)
      } else {
        object_contact_region(nose, obj_left, "radial", radius = contact_distance)
      }
      region_right <- if (box_side == "right") {
        object_contact_region(nose, obj_right, "box",
                              width = object_box_width, height = object_box_height)
      } else {
        object_contact_region(nose, obj_right, "radial", radius = contact_distance)
      }
    }
  } else {
    # One detector, both objects. Contact no longer depends on the novel-location
    # metadata, so a missing metadata row no longer discards all contact data.
    region_left <- object_contact_region(
      nose, obj_left, contact_geometry,
      radius = contact_distance, width = object_box_width, height = object_box_height
    )
    region_right <- object_contact_region(
      nose, obj_right, contact_geometry,
      radius = contact_distance, width = object_box_width, height = object_box_height
    )
  }

  contact_left <- region_left & body_left > body_exclusion_distance & oriented_left
  contact_right <- region_right & body_right > body_exclusion_distance & oriented_right

  proximity_left <- distance_left > proximity_range[1] & distance_left <= proximity_range[2]
  proximity_right <- distance_right > proximity_range[1] & distance_right <= proximity_range[2]
  proximity_angle_left <- proximity_left &
    abs(angle_left) >= proximity_angle[1] & abs(angle_left) <= proximity_angle[2]
  proximity_angle_right <- proximity_right &
    abs(angle_right) >= proximity_angle[1] & abs(angle_right) <= proximity_angle[2]

  segment <- function(event) {
    segment_events(event, fps, valid = valid, min_bout_s = min_bout_s, max_gap_s = max_gap_s)
  }
  left <- segment(contact_left)
  right <- segment(contact_right)
  proximity_left_seg <- segment(proximity_left)
  proximity_right_seg <- segment(proximity_right)
  proximity_angle_left_seg <- segment(proximity_angle_left)
  proximity_angle_right_seg <- segment(proximity_angle_right)

  # Preserve the repository's legacy assignment: metadata R maps the left side to
  # novel and metadata L maps the right side to novel. Its biological meaning must
  # be confirmed before changing it.
  novel_is_left <- !is.na(location) && location == "R"
  mapped <- function(left_value, right_value) {
    if (is.na(location)) return(c(novel = NA_real_, familiar = NA_real_))
    if (novel_is_left) c(novel = left_value, familiar = right_value) else
      c(novel = right_value, familiar = left_value)
  }
  contact_mapped <- mapped(left$duration_s, right$duration_s)
  proximity_mapped <- mapped(proximity_left_seg$duration_s, proximity_right_seg$duration_s)
  latency_values <- c(left$latency_s, right$latency_s)
  latency_values <- latency_values[is.finite(latency_values)]
  first_contact_latency <- if (length(latency_values) == 0) NA_real_ else min(latency_values)

  summary <- data.frame(
    contactLeft = left$duration_s,
    contactRight = right$duration_s,
    contactNov = unname(contact_mapped["novel"]),
    contactFam = unname(contact_mapped["familiar"]),
    proxLeft = proximity_left_seg$duration_s,
    proxRight = proximity_right_seg$duration_s,
    proxNov = unname(proximity_mapped["novel"]),
    proxFam = unname(proximity_mapped["familiar"]),
    proxLeftAngle = proximity_angle_left_seg$duration_s,
    proxRightAngle = proximity_angle_right_seg$duration_s,
    latency = first_contact_latency,
    latencyLeft = left$latency_s,
    latencyRight = right$latency_s,
    frequencyL = left$n_bouts,
    frequencyR = right$n_bouts,
    meanBoutLeft = left$mean_bout_s,
    meanBoutRight = right$mean_bout_s,
    entriesLeft = count_entries(contact_left, valid = valid),
    entriesRight = count_entries(contact_right, valid = valid),
    totalTime = get_tracking_duration(tracking),
    validTime = left$valid_time_s,
    validFraction = left$valid_time_s / get_tracking_duration(tracking),
    contactGeometry = contact_geometry,
    novelLoc = location,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    events = data.frame(
      contact_left = contact_left,
      contact_right = contact_right,
      proximity_left = proximity_left,
      proximity_right = proximity_right,
      valid = valid
    ),
    angles = data.frame(left = angle_left, right = angle_right),
    distances = data.frame(left = distance_left, right = distance_right),
    segmentation = list(left = left, right = right)
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
