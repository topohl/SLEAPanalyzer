# Generic behavioral event metrics derived from canonical bout tables.

event_duration <- function(event, fps, min_frames = 1L) {
  filtered <- filter_short_events(event, min_frames)
  frames_to_seconds(sum(filtered), fps)
}

event_frequency <- function(event, min_frames = 1L) {
  filtered <- filter_short_events(event, min_frames)
  if (length(filtered) == 0) return(0L)
  as.integer(sum(filtered & !c(FALSE, head(filtered, -1L))))
}

event_latency <- function(event, fps, min_frames = 1L) {
  bouts <- event_bout_table(event, fps, min_frames)
  if (nrow(bouts) == 0) return(NA_real_)
  bouts$latency_seconds[1]
}

summarize_event_metrics <- function(event, fps, min_frames = 1L) {
  filtered <- filter_short_events(event, min_frames)
  bouts <- event_bout_table(filtered, fps)
  list(
    duration_s = frames_to_seconds(sum(filtered), fps),
    bouts = as.integer(nrow(bouts)),
    latency_s = if (nrow(bouts) == 0) NA_real_ else bouts$latency_seconds[1],
    mean_bout_s = if (nrow(bouts) == 0) NA_real_ else mean(bouts$duration_seconds),
    max_bout_s = if (nrow(bouts) == 0) NA_real_ else max(bouts$duration_seconds),
    interbout_intervals_s = event_interbout_intervals(filtered, fps)
  )
}
