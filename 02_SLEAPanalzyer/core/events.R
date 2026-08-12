# Canonical logical-event handling and bout-table construction.

normalize_event_vector <- function(event) {
  event <- as.logical(event)
  event[is.na(event)] <- FALSE
  event
}

filter_short_events <- function(event, min_frames = 1L) {
  event <- normalize_event_vector(event)
  validate_positive_integer(min_frames, "min_frames")
  if (length(event) == 0 || min_frames == 1L) return(event)
  runs <- rle(event)
  runs$values <- runs$values & runs$lengths >= min_frames
  inverse.rle(runs)
}

empty_event_bout_table <- function() {
  data.frame(
    start_frame = integer(),
    end_frame = integer(),
    duration_frames = integer(),
    duration_seconds = numeric(),
    latency_seconds = numeric(),
    bout_number = integer()
  )
}

event_bout_table <- function(event, fps, min_frames = 1L) {
  validate_fps(fps)
  event <- filter_short_events(event, min_frames)
  if (!any(event)) return(empty_event_bout_table())
  starts <- which(event & !c(FALSE, head(event, -1L)))
  ends <- which(event & !c(tail(event, -1L), FALSE))
  durations <- ends - starts + 1L
  data.frame(
    start_frame = as.integer(starts),
    end_frame = as.integer(ends),
    duration_frames = as.integer(durations),
    duration_seconds = frames_to_seconds(durations, fps),
    latency_seconds = frames_to_seconds(starts - 1L, fps),
    bout_number = seq_along(starts)
  )
}

event_interbout_intervals <- function(event, fps, min_frames = 1L) {
  bouts <- event_bout_table(event, fps, min_frames)
  if (nrow(bouts) < 2) return(numeric())
  frames_to_seconds(
    bouts$start_frame[-1L] - bouts$end_frame[-nrow(bouts)] - 1L,
    fps
  )
}
