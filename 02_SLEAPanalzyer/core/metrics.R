# Generic behavioral event metrics derived from canonical bout tables.
#
# These are thin wrappers over segment_events(); prefer segment_events() plus
# event_summary_row() for new code because they also report valid analyzed
# time and percent-of-valid-time.

event_duration <- function(event, fps, min_frames = 1L, valid = NULL) {
  filtered <- filter_short_events(normalize_event_vector(event, valid), min_frames)
  frames_to_seconds(sum(filtered), fps)
}

event_frequency <- function(event, min_frames = 1L, valid = NULL) {
  filtered <- filter_short_events(normalize_event_vector(event, valid), min_frames)
  count_entries(filtered)
}

event_latency <- function(event, fps, min_frames = 1L, valid = NULL) {
  bouts <- event_bout_table(event, fps, min_frames, valid = valid)
  if (nrow(bouts) == 0) return(NA_real_)
  bouts$latency_seconds[1]
}

summarize_event_metrics <- function(event, fps, min_frames = 1L, valid = NULL,
                                    max_gap_frames = 0L) {
  validate_fps(fps)
  validate_positive_integer(min_frames, "min_frames")
  segmentation <- segment_events(
    event, fps,
    valid = valid,
    min_bout_frames = min_frames,
    max_gap_frames = max_gap_frames
  )
  list(
    duration_s = segmentation$duration_s,
    bouts = segmentation$n_bouts,
    latency_s = segmentation$latency_s,
    mean_bout_s = segmentation$mean_bout_s,
    median_bout_s = segmentation$median_bout_s,
    max_bout_s = segmentation$max_bout_s,
    valid_time_s = segmentation$valid_time_s,
    percent_valid_time = segmentation$percent_valid_time,
    interbout_intervals_s = segmentation$interbout_intervals_s
  )
}
