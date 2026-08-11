testthat::test_that("event summaries use frame time, bouts, and zero-based latency", {
  event <- c(rep(FALSE, 60), rep(TRUE, 30), rep(FALSE, 60))
  summary <- event_summary(event, fps = 30)
  testthat::expect_equal(summary$duration_s, 1)
  testthat::expect_equal(summary$bouts, 1L)
  testthat::expect_equal(summary$latency_s, 2)

  absent <- event_summary(rep(FALSE, 150), fps = 30)
  testthat::expect_true(is.na(absent$latency_s))
  testthat::expect_false(is.infinite(absent$latency_s))

  two <- c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE)
  testthat::expect_equal(event_summary(two, 2)$bouts, 2L)
  testthat::expect_equal(event_interbout_intervals_s(two, 2), 1)
})

testthat::test_that("short bouts and NA event values are handled deterministically", {
  event <- c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
  filtered <- suppress_short_event_bouts(event, min_frames = 2)
  testthat::expect_equal(filtered, c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))

  with_na <- c(FALSE, TRUE, NA, TRUE, FALSE)
  testthat::expect_equal(event_summary(with_na, fps = 1)$bouts, 2L)
  testthat::expect_equal(event_summary(with_na, fps = 1)$duration_s, 2)
})
