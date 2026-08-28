test_that("segment_events handles an event at the first frame", {
  event <- c(TRUE, TRUE, FALSE, FALSE)
  out <- segment_events(event, fps = 10)
  expect_equal(out$n_bouts, 1L)
  expect_equal(out$latency_s, 0)
  expect_equal(out$duration_s, 0.2)
  expect_equal(out$bouts$start_frame, 1L)
})

test_that("segment_events handles an event at the last frame", {
  event <- c(FALSE, FALSE, TRUE, TRUE)
  out <- segment_events(event, fps = 10)
  expect_equal(out$n_bouts, 1L)
  expect_equal(out$bouts$end_frame, 4L)
  expect_equal(out$latency_s, 0.2)
  expect_equal(out$duration_s, 0.2)
})

test_that("segment_events handles a one-frame event", {
  out <- segment_events(c(FALSE, TRUE, FALSE), fps = 25)
  expect_equal(out$n_bouts, 1L)
  expect_equal(out$duration_s, 1 / 25)
  expect_equal(out$mean_bout_s, 1 / 25)
  expect_equal(out$median_bout_s, 1 / 25)
})

test_that("no event yields censored NA latency, never Inf", {
  out <- segment_events(rep(FALSE, 50), fps = 30)
  expect_equal(out$n_bouts, 0L)
  expect_true(is.na(out$latency_s))
  expect_false(is.infinite(out$latency_s))
  expect_equal(out$duration_s, 0)
  expect_true(is.na(out$mean_bout_s))
  expect_equal(out$percent_valid_time, 0)
})

test_that("invalid frames are excluded from behavior and from analyzed time", {
  event <- rep(TRUE, 10)
  valid <- c(rep(TRUE, 6), rep(FALSE, 4))
  out <- segment_events(event, fps = 10, valid = valid)
  expect_equal(out$valid_frames, 6L)
  expect_equal(out$valid_time_s, 0.6)
  expect_equal(out$total_time_s, 1.0)
  expect_equal(out$duration_s, 0.6)
  expect_equal(out$percent_valid_time, 100)
})

test_that("an invalid gap inside an event splits the bout unless bridged", {
  event <- rep(TRUE, 11)
  valid <- c(rep(TRUE, 5), FALSE, rep(TRUE, 5))

  split <- segment_events(event, fps = 10, valid = valid)
  expect_equal(split$n_bouts, 2L)
  expect_equal(split$duration_s, 1.0)

  bridged <- segment_events(event, fps = 10, valid = valid, max_gap_s = 0.2)
  expect_equal(bridged$n_bouts, 1L)
  # Bridging joins the bouts but must not count the unobserved frame.
  expect_equal(bridged$duration_s, 1.0)
  expect_equal(bridged$valid_time_s, 1.0)
})

test_that("two bouts separated by a short gap merge only when allowed", {
  event <- c(rep(TRUE, 5), FALSE, FALSE, rep(TRUE, 5))
  expect_equal(segment_events(event, fps = 10)$n_bouts, 2L)
  expect_equal(segment_events(event, fps = 10, max_gap_s = 0.1)$n_bouts, 2L)
  expect_equal(segment_events(event, fps = 10, max_gap_s = 0.2)$n_bouts, 1L)
})

test_that("gap bridging never extends past the ends of the recording", {
  event <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  out <- segment_events(event, fps = 10, max_gap_s = 10)
  expect_equal(out$n_bouts, 1L)
  expect_equal(out$duration_s, 0.1)
})

test_that("all-invalid input reports no analyzed time and no behavior", {
  out <- segment_events(rep(TRUE, 20), fps = 30, valid = rep(FALSE, 20))
  expect_equal(out$valid_time_s, 0)
  expect_equal(out$duration_s, 0)
  expect_equal(out$n_bouts, 0L)
  expect_true(is.na(out$percent_valid_time))
  expect_true(is.na(out$latency_s))
})

test_that("NA in the event vector is treated as absence only when valid", {
  event <- c(TRUE, NA, TRUE)
  # Without a mask, NA collapses to FALSE (legacy behavior).
  expect_equal(segment_events(event, fps = 10)$n_bouts, 2L)
  # Deriving validity from the NA pattern marks the frame unobserved instead.
  out <- segment_events(event, fps = 10, valid = validity_from_event(event))
  expect_equal(out$valid_frames, 2L)
  expect_equal(out$valid_time_s, 0.2)
  expect_equal(out$duration_s, 0.2)
  expect_equal(out$percent_valid_time, 100)
})

test_that("minimum bout duration is expressed in seconds and scales with fps", {
  event <- c(rep(TRUE, 3), rep(FALSE, 5), rep(TRUE, 10))
  at_30 <- segment_events(event, fps = 30, min_bout_s = 0.2)
  expect_equal(at_30$n_bouts, 1L)
  at_10 <- segment_events(event, fps = 10, min_bout_s = 0.2)
  expect_equal(at_10$n_bouts, 2L)
})

test_that("minimum bout thresholds survive binary floating point", {
  # 0.1 * 30 is 3.0000000000000004; a naive ceiling would require four frames.
  event <- c(rep(FALSE, 2), rep(TRUE, 3), rep(FALSE, 2))
  expect_equal(segment_events(event, fps = 30, min_bout_s = 0.1)$n_bouts, 1L)
})

test_that("short bouts are filtered after gaps are bridged", {
  # Two two-frame bouts split by one frame: bridging first keeps them.
  event <- c(TRUE, TRUE, FALSE, TRUE, TRUE)
  unbridged <- segment_events(event, fps = 10, min_bout_s = 0.3)
  expect_equal(unbridged$n_bouts, 0L)
  bridged <- segment_events(event, fps = 10, min_bout_s = 0.3, max_gap_s = 0.1)
  expect_equal(bridged$n_bouts, 1L)
  expect_equal(bridged$duration_s, 0.5)
})

test_that("interbout intervals are reported in seconds", {
  event <- c(TRUE, FALSE, FALSE, FALSE, TRUE)
  out <- segment_events(event, fps = 10)
  expect_equal(out$interbout_intervals_s, 0.3)
})

test_that("bout summary statistics are consistent", {
  event <- c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
  out <- segment_events(event, fps = 10)
  expect_equal(out$n_bouts, 3L)
  expect_equal(out$duration_s, 0.6)
  expect_equal(out$max_bout_s, 0.3)
  expect_equal(out$median_bout_s, 0.2)
  expect_equal(out$mean_bout_s, 0.2)
})

test_that("count_entries counts onsets and count_state_changes counts both", {
  event <- c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
  expect_equal(count_entries(event), 2L)
  expect_equal(count_state_changes(event), 4L)

  # A run already in progress on frame one has no observed onset.
  started_inside <- c(TRUE, TRUE, FALSE)
  expect_equal(count_entries(started_inside), 0L)
  expect_equal(count_state_changes(started_inside), 1L)
})

test_that("a dropout inside a visit does not fabricate an entry", {
  event <- c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)
  expect_equal(count_entries(event), 1L)
  # Frame 4 is unobserved; the visit is still a single entry.
  expect_equal(count_entries(event, valid = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)), 1L)
})

test_that("count_initial exposes a visit already in progress at frame one", {
  event <- c(TRUE, TRUE, FALSE, TRUE)
  expect_equal(count_entries(event), 1L)
  expect_equal(count_entries(event, count_initial = TRUE), 2L)
})

test_that("event_summary_row flattens a segmentation", {
  out <- segment_events(c(TRUE, TRUE, FALSE), fps = 10)
  row <- event_summary_row(out, "center")
  expect_true(all(c(
    "center_duration_s", "center_percent_valid_time", "center_bouts",
    "center_latency_s", "center_valid_time_s"
  ) %in% names(row)))
  expect_equal(row$center_duration_s, 0.2)
  expect_equal(nrow(row), 1L)
})

test_that("validity masks reject length mismatches", {
  expect_error(
    segment_events(rep(TRUE, 5), fps = 10, valid = rep(TRUE, 4)),
    "validity mask has length"
  )
})

test_that("segment_events rejects invalid parameters", {
  expect_error(segment_events(rep(TRUE, 5), fps = 0), "fps")
  expect_error(segment_events(rep(TRUE, 5), fps = 10, min_bout_s = -1), "min_bout_s")
  expect_error(segment_events(rep(TRUE, 5), fps = 10, max_gap_s = -1), "max_gap_s")
})

test_that("summarize_event_metrics preserves frame-based minimum bouts", {
  event <- c(rep(TRUE, 3), FALSE, rep(TRUE, 5))
  out <- summarize_event_metrics(event, fps = 30, min_frames = 4L)
  expect_equal(out$bouts, 1L)
  expect_equal(out$duration_s, 5 / 30)
})
