test_that("short interior gaps are interpolated linearly and labelled", {
  x <- c(0, 1, NA, NA, 4, 5)
  out <- bounded_linear_interpolation(x, max_gap_frames = 2L)
  expect_equal(out$values, c(0, 1, 2, 3, 4, 5))
  expect_equal(
    as.character(out$status),
    c("observed", "observed", "interpolated", "interpolated", "observed", "observed")
  )
})

test_that("interior gaps longer than the limit stay missing", {
  x <- c(0, 1, NA, NA, NA, 5)
  out <- bounded_linear_interpolation(x, max_gap_frames = 2L)
  expect_true(all(is.na(out$values[3:5])))
  expect_equal(as.character(out$status[3:5]), rep("invalid", 3))
})

test_that("leading and trailing gaps are never filled", {
  x <- c(NA, NA, 2, 3, NA, NA)
  out <- bounded_linear_interpolation(x, max_gap_frames = 100L)
  expect_true(all(is.na(out$values[c(1, 2, 5, 6)])))
  expect_equal(
    as.character(out$status),
    c("invalid", "invalid", "observed", "observed", "invalid", "invalid")
  )
})

test_that("a zero gap limit disables interpolation entirely", {
  x <- c(0, NA, 2)
  out <- bounded_linear_interpolation(x, max_gap_frames = 0L)
  expect_true(is.na(out$values[2]))
  expect_equal(as.character(out$status[2]), "invalid")
})

test_that("an all-missing vector stays missing", {
  out <- bounded_linear_interpolation(rep(NA_real_, 5), max_gap_frames = 10L)
  expect_true(all(is.na(out$values)))
  expect_equal(as.character(out$status), rep("invalid", 5))
})

test_that("a supplied validity mask overrides finite but rejected values", {
  x <- c(0, 99, 2)
  out <- bounded_linear_interpolation(x, valid = c(TRUE, FALSE, TRUE), max_gap_frames = 1L)
  expect_equal(out$values, c(0, 1, 2))
  expect_equal(as.character(out$status[2]), "interpolated")
})

test_that("interpolation is exactly linear across a multi-frame gap", {
  x <- c(10, rep(NA_real_, 4), 20)
  out <- bounded_linear_interpolation(x, max_gap_frames = 4L)
  expect_equal(out$values, c(10, 12, 14, 16, 18, 20))
})

test_that("interpolate_tracking uses seconds and records status", {
  n <- 10
  points <- list(
    nose = cbind(x = as.numeric(1:n), y = as.numeric(1:n))
  )
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  tracking$data$nose$x[5] <- NA_real_
  tracking$data$nose$y[5] <- NA_real_
  tracking$data$nose$x[8:9] <- NA_real_
  tracking$data$nose$y[8:9] <- NA_real_

  # 0.1 s at 10 fps is a one-frame limit: fills frame 5, leaves 8-9 missing.
  filled <- interpolate_tracking(tracking, "nose", max_gap_s = 0.1)
  expect_equal(filled$data$nose$x[5], 5)
  expect_true(all(is.na(filled$data$nose$x[8:9])))
  expect_equal(as.character(filled$data$nose$status[5]), "interpolated")
  expect_equal(as.character(filled$data$nose$status[8:9]), c("invalid", "invalid"))
  expect_equal(filled$interpolation$max_gap_frames, 1L)
})

test_that("x and y always share one status decision", {
  n <- 6
  points <- list(nose = cbind(x = as.numeric(1:n), y = as.numeric(1:n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  # Only y is missing; the frame must still be treated as unobserved for both.
  tracking$data$nose$y[3] <- NA_real_

  filled <- interpolate_tracking(tracking, "nose", max_gap_s = 0)
  expect_true(is.na(filled$data$nose$x[3]))
  expect_true(is.na(filled$data$nose$y[3]))
  expect_equal(as.character(filled$data$nose$status[3]), "invalid")
})

test_that("a likelihood cutoff rejects low-confidence frames", {
  n <- 6
  points <- list(nose = cbind(x = as.numeric(1:n), y = as.numeric(1:n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  tracking$data$nose$likelihood[4] <- 0.1

  filled <- interpolate_tracking(tracking, "nose", max_gap_s = 0.1, likelihood_cutoff = 0.9)
  expect_equal(as.character(filled$data$nose$status[4]), "interpolated")
  expect_equal(filled$data$nose$x[4], 4)

  strict <- interpolate_tracking(tracking, "nose", max_gap_s = 0, likelihood_cutoff = 0.9)
  expect_true(is.na(strict$data$nose$x[4]))
})

test_that("landmark_validity treats interpolated frames as invalid by default", {
  n <- 6
  points <- list(nose = cbind(x = as.numeric(1:n), y = as.numeric(1:n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  tracking$data$nose$x[3] <- NA_real_
  tracking$data$nose$y[3] <- NA_real_
  filled <- interpolate_tracking(tracking, "nose", max_gap_s = 0.1)

  expect_equal(sum(landmark_validity(filled, "nose")), 5)
  expect_equal(
    sum(landmark_validity(filled, "nose", interpolated_is_valid = TRUE)), 6
  )
})

test_that("interpolation_report summarises provenance per landmark", {
  n <- 10
  points <- list(nose = cbind(x = as.numeric(1:n), y = as.numeric(1:n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  tracking$data$nose$x[4] <- NA_real_
  tracking$data$nose$y[4] <- NA_real_
  tracking$data$nose$x[7:9] <- NA_real_
  tracking$data$nose$y[7:9] <- NA_real_
  filled <- interpolate_tracking(tracking, "nose", max_gap_s = 0.1)

  report <- interpolation_report(filled, "nose")
  expect_equal(report$observed_frames, 6L)
  expect_equal(report$interpolated_frames, 1L)
  expect_equal(report$invalid_frames, 3L)
  expect_equal(report$longest_invalid_gap_frames, 3L)
  expect_equal(report$longest_invalid_gap_s, 0.3)
})
