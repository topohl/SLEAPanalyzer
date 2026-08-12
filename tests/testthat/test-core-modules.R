testthat::test_that("TrackingData accessors validate shape, landmarks, fps, and alignment", {
  tracking <- make_tracking(
    frames = 10:12, fps = 20,
    points = list(nose = cbind(x = 1:3, y = 4:6))
  )
  testthat::expect_true(validate_tracking_data(tracking))
  testthat::expect_equal(get_tracking_frames(tracking), 10:12)
  testthat::expect_equal(get_tracking_fps(tracking), 20)
  testthat::expect_equal(get_tracking_duration(tracking), 3 / 20)
  testthat::expect_true(has_landmarks(tracking, "nose"))
  testthat::expect_false(has_landmarks(tracking, "tail"))
  testthat::expect_error(
    validate_tracking_data(tracking, required_landmarks = "tail"),
    "missing landmark"
  )

  invalid <- tracking
  invalid$object.type <- "not-tracking"
  testthat::expect_error(validate_tracking_data(invalid), "not a canonical")
  invalid <- tracking
  invalid$fps <- 0
  testthat::expect_error(validate_tracking_data(invalid), "fps must be positive")
  invalid <- tracking
  invalid$data$nose$frame <- 11:13
  testthat::expect_error(validate_tracking_data(invalid), "not aligned")
})

testthat::test_that("coordinate and threshold units are explicit and conversions preserve units", {
  testthat::expect_equal(canonical_coordinate_unit("pixels"), "px")
  testthat::expect_equal(canonical_coordinate_unit("cm"), "cm")
  testthat::expect_equal(validate_coordinate_unit("pixels"), "pixels")
  testthat::expect_equal(validate_coordinate_unit("cm"), "cm")
  testthat::expect_equal(frames_to_seconds(c(0, 15, 30), 30), c(0, 0.5, 1))
  testthat::expect_equal(distance_to_speed(c(0, 2), 10), c(0, 20))
  testthat::expect_equal(speed_unit("pixels"), "px/s")
  testthat::expect_error(validate_threshold_unit("cm", "pixel"), "does not match")

  tracking <- make_tracking(distance_unit = "pixel")
  testthat::expect_error(
    compute_nor_metrics(tracking, "R", fps = 30),
    "missing landmark"
  )
})

testthat::test_that("geometry primitives return known distances, angles, and membership", {
  testthat::expect_equal(euclidean_distance(0, 0, 3, 4), 5)
  testthat::expect_equal(vector_angle_degrees(1, 0, 0, 1), 90)
  testthat::expect_equal(vector_angle_degrees(1, 0, -1, 0), 180)
  testthat::expect_true(is.na(vector_angle_degrees(0, 0, 1, 0)))

  polygon <- data.frame(x = c(0, 2, 2, 0), y = c(0, 0, 2, 2))
  testthat::expect_equal(polygon_area(polygon), 4)
  testthat::expect_equal(
    points_in_polygon(c(1, 3, 0), c(1, 1, 1), polygon),
    c(TRUE, FALSE, TRUE)
  )
  testthat::expect_equal(
    points_in_axis_aligned_box(c(0, 2), c(0, 0), c(0, 0), c(0, 0), 2, 2),
    c(TRUE, FALSE)
  )
})

testthat::test_that("canonical bout tables cover boundaries, gaps, NA, and minimum durations", {
  single <- event_bout_table(c(FALSE, TRUE, TRUE, FALSE), fps = 2)
  testthat::expect_equal(
    single,
    data.frame(
      start_frame = 2L, end_frame = 3L, duration_frames = 2L,
      duration_seconds = 1, latency_seconds = 0.5, bout_number = 1L
    )
  )

  two <- event_bout_table(c(TRUE, FALSE, FALSE, TRUE), fps = 2)
  testthat::expect_equal(two$start_frame, c(1L, 4L))
  testthat::expect_equal(two$end_frame, c(1L, 4L))
  testthat::expect_equal(two$latency_seconds, c(0, 1.5))
  testthat::expect_equal(event_interbout_intervals(c(TRUE, FALSE, FALSE, TRUE), 2), 1)

  absent <- event_bout_table(c(FALSE, NA, FALSE), fps = 30)
  testthat::expect_equal(nrow(absent), 0)
  testthat::expect_named(
    absent,
    c("start_frame", "end_frame", "duration_frames", "duration_seconds",
      "latency_seconds", "bout_number")
  )

  filtered <- event_bout_table(
    c(TRUE, FALSE, TRUE, TRUE, TRUE), fps = 1, min_frames = 2
  )
  testthat::expect_equal(filtered$start_frame, 3L)
  testthat::expect_equal(filtered$end_frame, 5L)
})

testthat::test_that("QC helpers report missingness and likelihood without changing data", {
  tracking <- make_tracking(
    frames = 0:2,
    points = list(nose = cbind(x = c(0, NA, 2), y = c(0, 1, 2)))
  )
  before <- tracking
  report <- tracking_qc_report(tracking, likelihood_cutoff = 0.8)
  testthat::expect_equal(report$coordinates$valid_frames, 2)
  testthat::expect_equal(report$coordinates$missing_or_invalid_frames, 1)
  testthat::expect_equal(report$likelihood$below_cutoff, 0L)
  testthat::expect_identical(tracking, before)

  degenerate <- data.frame(x = 0:3, y = 0:3)
  testthat::expect_false(calibration_qc_report(degenerate, 10, 10)$valid)
})

testthat::test_that("NOR validates timing and coordinate units before applying thresholds", {
  n <- 3
  points <- list(
    bodycentre = cbind(x = rep(-1, n), y = rep(0, n)),
    nose = cbind(x = rep(0, n), y = rep(0, n)),
    objL = cbind(x = rep(2, n), y = rep(0, n)),
    objR = cbind(x = rep(-3, n), y = rep(0, n))
  )
  pixels <- make_tracking(frames = 0:(n - 1), points = points, distance_unit = "pixel")
  testthat::expect_error(compute_nor_metrics(pixels, "R", 30), "does not match")

  metric <- make_tracking(frames = 0:(n - 1), points = points, distance_unit = "cm")
  testthat::expect_error(compute_nor_metrics(metric, "R", 10), "does not match")
})
