# Canonical logical-event handling and bout-table construction.
#
# Every assay must segment behavioral events through this module. The engine is
# validity-aware: frames that were never observed are neither counted as
# behavior nor as analyzed time. See docs/event_engine.md.

# ---------------------------------------------------------------------------
# Event / validity vectors
# ---------------------------------------------------------------------------

#' Coerce an event vector to a strict logical vector.
#'
#' NA means "unknown", not "no behavior". Callers that have a validity mask
#' must pass it so unknown frames are recorded as invalid rather than silently
#' scored as absence of the behavior.
#'
#' @param event a logical (or coercible) vector
#' @param valid optional logical validity mask of the same length
#' @return a logical vector with no NA values
normalize_event_vector <- function(event, valid = NULL) {
  event <- as.logical(event)
  if (!is.null(valid)) {
    valid <- normalize_validity_mask(valid, length(event))
    event[!valid] <- FALSE
  }
  event[is.na(event)] <- FALSE
  event
}

#' Build a strict logical validity mask.
#'
#' NA in a validity mask is treated as invalid, which is the conservative
#' direction: an unknown frame is never counted as analyzed time.
#'
#' @param valid a logical vector, or NULL for "all frames valid"
#' @param n the required length
#' @return a logical vector of length n with no NA values
normalize_validity_mask <- function(valid, n) {
  validate_nonnegative_integer(n, "n")
  if (is.null(valid)) return(rep(TRUE, n))
  valid <- as.logical(valid)
  if (length(valid) != n) {
    stop("validity mask has length ", length(valid), " but ", n, " frames were expected")
  }
  valid[is.na(valid)] <- FALSE
  valid
}

#' Derive a validity mask from the NA pattern of an event vector.
#'
#' Useful for legacy call sites that only have the raw comparison result.
#'
#' @param event a logical (or coercible) vector
#' @return a logical vector that is FALSE wherever event is NA
validity_from_event <- function(event) {
  !is.na(as.logical(event))
}

# ---------------------------------------------------------------------------
# Run-length shaping
# ---------------------------------------------------------------------------

#' Drop event runs shorter than min_frames.
filter_short_events <- function(event, min_frames = 1L) {
  event <- normalize_event_vector(event)
  validate_positive_integer(min_frames, "min_frames")
  if (length(event) == 0 || min_frames == 1L) return(event)
  runs <- rle(event)
  runs$values <- runs$values & runs$lengths >= min_frames
  inverse.rle(runs)
}

#' Bridge non-event runs no longer than max_gap_frames.
#'
#' Only interior gaps are bridged; leading and trailing non-event runs are
#' never converted into behavior.
#'
#' @param event a logical vector
#' @param max_gap_frames maximum bridged gap in frames (0 disables merging)
#' @return a logical vector
merge_event_gaps <- function(event, max_gap_frames = 0L) {
  event <- normalize_event_vector(event)
  validate_nonnegative_integer(max_gap_frames, "max_gap_frames")
  if (length(event) == 0 || max_gap_frames == 0L || !any(event)) return(event)
  runs <- rle(event)
  interior <- seq_along(runs$values)
  interior <- interior[interior != 1L & interior != length(runs$values)]
  bridge <- interior[!runs$values[interior] & runs$lengths[interior] <= max_gap_frames]
  runs$values[bridge] <- TRUE
  inverse.rle(runs)
}

# ---------------------------------------------------------------------------
# Bout tables
# ---------------------------------------------------------------------------

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

#' Tabulate the bouts of a logical event vector.
#'
#' start_frame / end_frame are 1-based positions in event, not raw video frame
#' identifiers. latency_seconds is measured from the first analyzed frame.
#'
#' @param event a logical vector
#' @param fps frames per second
#' @param min_frames minimum bout length in frames, applied after merging
#' @param valid optional validity mask; invalid frames cannot be event frames
#' @param max_gap_frames optional interior gap bridged before filtering
event_bout_table <- function(event, fps, min_frames = 1L, valid = NULL,
                             max_gap_frames = 0L) {
  validate_fps(fps)
  event <- normalize_event_vector(event, valid)
  event <- merge_event_gaps(event, max_gap_frames)
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

#' Intervals between consecutive bouts, in seconds.
event_interbout_intervals <- function(event, fps, min_frames = 1L, valid = NULL,
                                      max_gap_frames = 0L) {
  bouts <- event_bout_table(event, fps, min_frames, valid, max_gap_frames)
  if (nrow(bouts) < 2) return(numeric())
  frames_to_seconds(
    bouts$start_frame[-1L] - bouts$end_frame[-nrow(bouts)] - 1L,
    fps
  )
}

# ---------------------------------------------------------------------------
# Canonical segmentation entry point
# ---------------------------------------------------------------------------

#' Segment a behavioral event using seconds-valued parameters.
#'
#' This is the canonical entry point for every assay. Durations are expressed
#' in seconds and converted to frames with fps, so the same configuration is
#' valid at any acquisition rate.
#'
#' Semantics:
#'   - Invalid frames are never event frames and never count as analyzed time.
#'   - max_gap_s bridges interior gaps before min_bout_s is applied, so a bout
#'     briefly interrupted by a dropout is not split into two short bouts that
#'     both fall below the minimum. Bridging never marks an invalid frame as
#'     behavior.
#'   - When no bout occurs, latency_s is NA_real_ (right-censored), never Inf
#'     and never the trial duration.
#'
#' @param event a logical vector of per-frame behavior
#' @param fps frames per second
#' @param valid optional per-frame validity mask
#' @param min_bout_s minimum bout duration in seconds
#' @param max_gap_s maximum bridged interior gap in seconds
#' @param min_bout_frames optional exact frame count overriding min_bout_s
#' @param max_gap_frames optional exact frame count overriding max_gap_s
#' @return a list with the bout table and the standard summary quantities
segment_events <- function(event, fps, valid = NULL, min_bout_s = 0,
                           max_gap_s = 0, min_bout_frames = NULL,
                           max_gap_frames = NULL) {
  validate_fps(fps)

  n <- length(event)
  valid <- normalize_validity_mask(valid, n)

  if (is.null(min_bout_frames)) {
    validate_scalar_number(min_bout_s, "min_bout_s", positive = TRUE, allow_zero = TRUE)
    # Round before ceiling: 0.1 * 30 is 3.0000000000000004 in binary floating
    # point, which would otherwise demand a four-frame minimum.
    min_frames <- max(1L, as.integer(ceiling(round(min_bout_s * fps, 9))))
  } else {
    min_frames <- validate_positive_integer(min_bout_frames, "min_bout_frames")
    min_bout_s <- min_frames / fps
  }
  if (is.null(max_gap_frames)) {
    validate_scalar_number(max_gap_s, "max_gap_s", positive = TRUE, allow_zero = TRUE)
    max_gap_frames <- as.integer(floor(round(max_gap_s * fps, 9)))
  } else {
    max_gap_frames <- validate_nonnegative_integer(max_gap_frames, "max_gap_frames")
    max_gap_s <- max_gap_frames / fps
  }

  # Two vectors are kept deliberately.
  #
  # `segmented` defines bout *identity*: a behavior briefly interrupted by a
  # tracking dropout is one episode, not two. `observed` defines what was
  # actually *seen*: it never includes a frame that was not observed, so
  # durations and percentages are never inflated by bridged gaps.
  segmented <- normalize_event_vector(event, valid)
  segmented <- merge_event_gaps(segmented, max_gap_frames)
  segmented <- filter_short_events(segmented, min_frames)
  observed <- segmented & valid

  bouts <- event_bout_table(segmented, fps)
  if (nrow(bouts) > 0) {
    bouts$observed_frames <- vapply(
      seq_len(nrow(bouts)),
      function(i) sum(observed[bouts$start_frame[i]:bouts$end_frame[i]]),
      integer(1)
    )
    bouts$observed_seconds <- frames_to_seconds(bouts$observed_frames, fps)
  } else {
    bouts$observed_frames <- integer()
    bouts$observed_seconds <- numeric()
  }

  valid_time_s <- frames_to_seconds(sum(valid), fps)
  duration_s <- frames_to_seconds(sum(observed), fps)

  list(
    event = segmented,
    observed_event = observed,
    valid = valid,
    bouts = bouts,
    fps = fps,
    min_bout_s = min_bout_s,
    max_gap_s = max_gap_s,
    total_frames = n,
    valid_frames = sum(valid),
    valid_time_s = valid_time_s,
    total_time_s = frames_to_seconds(n, fps),
    duration_s = duration_s,
    percent_valid_time = if (valid_time_s > 0) duration_s / valid_time_s * 100 else NA_real_,
    n_bouts = as.integer(nrow(bouts)),
    latency_s = if (nrow(bouts) == 0) NA_real_ else bouts$latency_seconds[1],
    mean_bout_s = if (nrow(bouts) == 0) NA_real_ else mean(bouts$observed_seconds),
    median_bout_s = if (nrow(bouts) == 0) NA_real_ else stats::median(bouts$observed_seconds),
    max_bout_s = if (nrow(bouts) == 0) NA_real_ else max(bouts$observed_seconds),
    interbout_intervals_s = if (nrow(bouts) < 2) numeric() else {
      frames_to_seconds(
        bouts$start_frame[-1L] - bouts$end_frame[-nrow(bouts)] - 1L, fps
      )
    }
  )
}

#' Flatten segment_events() into a one-row data frame.
#'
#' @param segmentation the result of segment_events()
#' @param prefix column-name prefix, for example "center"
event_summary_row <- function(segmentation, prefix) {
  if (length(prefix) != 1 || !is.character(prefix) || !nzchar(prefix)) {
    stop("prefix must be one non-empty character value")
  }
  out <- data.frame(
    duration_s = segmentation$duration_s,
    percent_valid_time = segmentation$percent_valid_time,
    bouts = segmentation$n_bouts,
    latency_s = segmentation$latency_s,
    mean_bout_s = segmentation$mean_bout_s,
    median_bout_s = segmentation$median_bout_s,
    max_bout_s = segmentation$max_bout_s,
    valid_time_s = segmentation$valid_time_s,
    stringsAsFactors = FALSE
  )
  names(out) <- paste(prefix, names(out), sep = "_")
  out
}

# ---------------------------------------------------------------------------
# Entries and transitions
# ---------------------------------------------------------------------------

#' Count zone entries (observed onsets only).
#'
#' Distinct from the legacy CalculateTransitions(), which counted onsets *and*
#' offsets and therefore returned roughly twice the entry count, with a bias
#' depending on whether the animal started or ended inside the zone.
#'
#' Two properties matter for correctness:
#'
#' 1. Unobserved frames carry no state. Entries are detected between
#'    consecutive *observed* frames, so a tracking dropout in the middle of a
#'    zone visit does not fabricate an extra entry.
#' 2. A visit already in progress on the first observed frame is not an entry,
#'    because its onset was never observed. Set `count_initial = TRUE` to score
#'    it as an entry instead (for example when the animal is deliberately
#'    released into the zone).
#'
#' Note that this differs from the bout count reported by `segment_events()`:
#' an episode in progress at the first frame is still an observed bout, it
#' simply has no observed onset.
#'
#' @param event a logical vector
#' @param valid optional validity mask
#' @param count_initial whether a visit in progress at the first observed frame
#'   counts as an entry
#' @return an integer entry count
count_entries <- function(event, valid = NULL, count_initial = FALSE) {
  event <- as.logical(event)
  valid <- normalize_validity_mask(valid, length(event))
  observed <- event[valid]
  observed[is.na(observed)] <- FALSE
  if (length(observed) == 0) return(0L)
  onsets <- observed & !c(FALSE, head(observed, -1L))
  if (!isTRUE(count_initial)) onsets[1] <- FALSE
  as.integer(sum(onsets))
}

#' Count state changes (onsets plus offsets) between consecutive observations.
#'
#' Retained so that summaries produced with the legacy `CalculateTransitions()`
#' definition can still be reproduced. Prefer `count_entries()` for new work:
#' this quantity is roughly twice the entry count and depends on the state at
#' the first and last frame.
count_state_changes <- function(event, valid = NULL) {
  event <- as.logical(event)
  valid <- normalize_validity_mask(valid, length(event))
  observed <- event[valid]
  observed[is.na(observed)] <- FALSE
  if (length(observed) < 2) return(0L)
  as.integer(sum(observed[-1L] != observed[-length(observed)]))
}
