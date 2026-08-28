# Canonical QC pipeline: one synthetic failure mode per test.

qc_tracking <- function(n = 30, fps = 10, ...) {
  points <- list(...)
  if (length(points) == 0) {
    points <- list(bodycentre = cbind(x = as.numeric(1:n), y = rep(0, n)))
  }
  make_tracking(frames = 0:(n - 1), fps = fps, points = points)
}

test_that("a clean recording reports full valid time and passes", {
  tracking <- qc_tracking()
  report <- tracking_qc_report(tracking)
  expect_equal(report$valid_fraction, 1)
  expect_equal(report$valid_time_s, 3)
  expect_equal(report$longest_invalid_gap_s, 0)
  expect_true(qc_flags(report)$pass)
})

test_that("missing coordinates reduce valid time and are counted", {
  tracking <- qc_tracking()
  tracking$data$bodycentre$x[5:14] <- NA_real_
  tracking$data$bodycentre$y[5:14] <- NA_real_

  report <- tracking_qc_report(tracking)
  expect_equal(report$valid_frames, 20L)
  expect_equal(report$valid_time_s, 2)
  expect_equal(report$valid_fraction, 2 / 3)
  expect_equal(report$longest_invalid_gap_frames, 10L)
  expect_equal(report$longest_invalid_gap_s, 1)
  expect_equal(report$coordinates$missing_or_invalid_frames, 10L)

  flags <- qc_flags(report, min_valid_fraction = 0.8)
  expect_false(flags$pass)
  expect_match(flags$reasons[1], "valid fraction")
})

test_that("low confidence is reported as a fraction", {
  tracking <- qc_tracking()
  tracking$data$bodycentre$likelihood[1:6] <- 0.1

  report <- tracking_qc_report(tracking, likelihood_cutoff = 0.9)
  expect_equal(report$likelihood$below_cutoff, 6L)
  expect_equal(report$likelihood$low_confidence_fraction, 0.2)
  # Low-confidence frames do not count as analyzed time.
  expect_equal(report$valid_frames, 24L)
})

test_that("the longest invalid gap is distinguished from total missingness", {
  tracking <- qc_tracking()
  # Same total missing, very different consequences.
  scattered <- tracking
  scattered$data$bodycentre$x[seq(2, 20, by = 2)] <- NA_real_
  scattered$data$bodycentre$y[seq(2, 20, by = 2)] <- NA_real_
  contiguous <- tracking
  contiguous$data$bodycentre$x[1:10] <- NA_real_
  contiguous$data$bodycentre$y[1:10] <- NA_real_

  a <- tracking_qc_report(scattered)
  b <- tracking_qc_report(contiguous)
  expect_equal(a$valid_frames, b$valid_frames)
  expect_equal(a$longest_invalid_gap_frames, 1L)
  expect_equal(b$longest_invalid_gap_frames, 10L)

  # Both have the same valid fraction, so relax that criterion to isolate the
  # gap-length criterion: scattered dropouts pass, one long gap does not.
  expect_true(qc_flags(a, min_valid_fraction = 0.5, max_longest_gap_s = 0.5)$pass)
  gap_flags <- qc_flags(b, min_valid_fraction = 0.5, max_longest_gap_s = 0.5)
  expect_false(gap_flags$pass)
  expect_length(gap_flags$reasons, 1)
  expect_match(gap_flags$reasons, "longest gap")
})

test_that("interpolated time is separated from observed time", {
  tracking <- qc_tracking()
  tracking$data$bodycentre$x[10] <- NA_real_
  tracking$data$bodycentre$y[10] <- NA_real_
  filled <- interpolate_tracking(tracking, "bodycentre", max_gap_s = 0.5)

  report <- tracking_qc_report(filled)
  expect_equal(report$observed_time_s, 2.9)
  expect_equal(report$interpolated_time_s, 0.1)
  expect_equal(report$valid_time_s, 3)
  expect_equal(report$provenance$interpolated_frames, 1L)

  expect_false(qc_flags(report, max_interpolated_fraction = 0.01)$pass)
  expect_true(qc_flags(report, max_interpolated_fraction = 0.5)$pass)
})

test_that("implausible displacement is flagged against a speed limit", {
  n <- 20
  x <- as.numeric(1:n)
  x[10] <- 5000
  tracking <- qc_tracking(n = n, bodycentre = cbind(x = x, y = rep(0, n)))

  # Normal motion is 1 unit per frame at 10 fps, so 10 units/s.
  report <- tracking_qc_report(tracking, max_speed = 100)
  expect_equal(report$displacement$implausible_intervals, 2L)
  expect_gt(report$displacement$max_speed, 1000)
  expect_equal(report$displacement$median_speed, 10)

  no_limit <- tracking_qc_report(tracking)
  expect_true(is.na(no_limit$displacement$implausible_intervals))
})

test_that("abnormal body length is detected", {
  n <- 20
  nose <- cbind(x = as.numeric(1:n) + 5, y = rep(0, n))
  tail <- cbind(x = as.numeric(1:n), y = rep(0, n))
  # One frame where the nose is labelled far away.
  nose[10, "x"] <- 100
  tracking <- qc_tracking(n = n, nose = nose, tailbase = tail)

  report <- skeleton_qc(tracking, "nose", "tailbase", tolerance = 0.5)
  expect_equal(report$median_length, 5)
  expect_equal(report$abnormal_frames, 1L)
  expect_equal(report$abnormal_fraction, 1 / 20)
})

test_that("an identity swap is distinguished from two fast animals", {
  n <- 10
  # Animals sitting still, then exchanging positions on frame 6.
  ax <- c(rep(0, 5), rep(100, 5))
  bx <- c(rep(100, 5), rep(0, 5))
  swapped <- qc_tracking(
    n = n,
    a = cbind(x = ax, y = rep(0, n)),
    b = cbind(x = bx, y = rep(0, n))
  )
  report <- identity_swap_qc(swapped, "a", "b", jump_threshold = 50)
  expect_equal(report$suspected_swaps, 1L)

  # Both animals sprinting in the same direction: large jumps, no swap,
  # because neither lands where the other was.
  moving <- qc_tracking(
    n = n,
    a = cbind(x = seq(0, by = 100, length.out = n), y = rep(0, n)),
    b = cbind(x = seq(500, by = 100, length.out = n), y = rep(0, n))
  )
  clean <- identity_swap_qc(moving, "a", "b", jump_threshold = 50)
  expect_equal(clean$suspected_swaps, 0L)
})

test_that("landmark frame-count consistency is reported", {
  tracking <- qc_tracking(
    n = 10,
    a = cbind(x = rep(0, 10), y = rep(0, 10)),
    b = cbind(x = rep(1, 10), y = rep(1, 10))
  )
  expect_true(frame_consistency_qc(tracking)$consistent)
  expect_equal(frame_consistency_qc(tracking)$landmarks, 2L)
})

test_that("arena violations appear in the combined report once calibrated", {
  n <- 10
  corners <- data.frame(x = c(0, 100, 100, 0), y = c(100, 100, 0, 0))
  calibration <- arena_calibration(corners, 50, 50)
  x <- rep(50, n)
  x[3] <- 100000
  tracking <- qc_tracking(n = n, bodycentre = cbind(x = x, y = rep(50, n)))
  rectified <- rectify_tracking(tracking, calibration)

  report <- tracking_qc_report(rectified)
  expect_equal(report$arena$outside_arena_frames, 1L)
  expect_equal(report$coordinate_unit, "cm")
})

test_that("qc_summary_row flattens a report for batch tables", {
  tracking <- qc_tracking()
  row <- qc_summary_row(tracking_qc_report(tracking))
  expect_equal(nrow(row), 1L)
  expect_true(all(c(
    "file", "fps", "valid_time_s", "valid_fraction", "longest_invalid_gap_s",
    "coordinate_unit"
  ) %in% names(row)))
  expect_equal(row$valid_time_s, 3)
})

test_that("qc_flags reports every failing reason without dropping data", {
  tracking <- qc_tracking()
  tracking$data$bodycentre$x[1:20] <- NA_real_
  tracking$data$bodycentre$y[1:20] <- NA_real_
  report <- tracking_qc_report(tracking)

  flags <- qc_flags(report, min_valid_fraction = 0.9, max_longest_gap_s = 0.5)
  expect_false(flags$pass)
  expect_length(flags$reasons, 2)
  # QC is report-only: the tracking object is untouched.
  expect_equal(report$frames, 30L)
})
