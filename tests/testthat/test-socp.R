# Synthetic Social Preference trajectories with known expected results.

socp_tracking <- function(n = 60, fps = 30, nose, bodycentre, socl, socr,
                          unit = "cm") {
  points <- square_points(n)
  rep_xy <- function(xy) {
    if (is.matrix(xy)) return(xy)
    cbind(x = rep(xy[1], n), y = rep(xy[2], n))
  }
  points$nose <- rep_xy(nose)
  points$bodycentre <- rep_xy(bodycentre)
  points$socl <- rep_xy(socl)
  points$socr <- rep_xy(socr)
  make_tracking(frames = 0:(n - 1), fps = fps, points = points, distance_unit = unit)
}

test_that("SocP rejects centimetre thresholds on pixel coordinates", {
  tracking <- socp_tracking(
    n = 10, fps = 10, nose = c(0, 0), bodycentre = c(-5, 0),
    socl = c(2, 0), socr = c(100, 0), unit = "pixel"
  )
  # The pre-v2 implementation applied the 6 cm contact threshold to pixel
  # coordinates without complaint.
  expect_error(
    compute_socp_metrics(tracking, "R", fps = 10, threshold_unit = "cm"),
    "does not match"
  )
  expect_silent(compute_socp_metrics(tracking, "R", fps = 10))
})

test_that("SocP rejects an fps that disagrees with the tracking data", {
  tracking <- socp_tracking(
    n = 10, fps = 10, nose = c(0, 0), bodycentre = c(-5, 0),
    socl = c(2, 0), socr = c(100, 0)
  )
  expect_error(compute_socp_metrics(tracking, "R", fps = 30), "fps does not match")
})

test_that("SocP requires the stimulus landmarks", {
  tracking <- socp_tracking(
    n = 10, fps = 10, nose = c(0, 0), bodycentre = c(-5, 0),
    socl = c(2, 0), socr = c(100, 0)
  )
  tracking$data$socr <- NULL
  expect_error(compute_socp_metrics(tracking, "R", fps = 10), "socr")
})

test_that("the same detector is applied to both chambers", {
  n <- 30
  # The nose sits exactly halfway between the two stimulus chambers.
  tracking <- socp_tracking(
    n = n, fps = 30, nose = c(0, 0), bodycentre = c(0, -5),
    socl = c(-5, 0), socr = c(5, 0)
  )
  result <- compute_socp_metrics(tracking, "R", fps = 30)
  expect_equal(result$summary$contactLeft, result$summary$contactRight)
  expect_equal(result$summary$contactNovel, result$summary$contactFamiliar)
})

test_that("missing frames are excluded from behavior and from analyzed time", {
  n <- 30
  fps <- 10
  tracking <- socp_tracking(
    n = n, fps = fps, nose = c(2, 0), bodycentre = c(-5, 0),
    socl = c(2, 0), socr = c(100, 0)
  )
  tracking$data$nose$x[11:20] <- NA_real_
  tracking$data$nose$y[11:20] <- NA_real_

  result <- compute_socp_metrics(tracking, "R", fps = fps)
  expect_equal(result$summary$totalTime, 3)
  expect_equal(result$summary$validTime, 2)
  expect_equal(result$summary$validFraction, 2 / 3)
  # Pre-v2, the ten untracked frames were scored as confident non-contact and
  # still counted toward the denominator.
  expect_equal(result$summary$contactLeft, 2)
})

test_that("novel and familiar follow the metadata, and missing metadata is safe", {
  n <- 20
  fps <- 10
  tracking <- socp_tracking(
    n = n, fps = fps, nose = c(2, 0), bodycentre = c(-5, 0),
    socl = c(2, 0), socr = c(100, 0)
  )

  right_novel <- compute_socp_metrics(tracking, "L", fps = fps)
  expect_equal(right_novel$summary$contactNovel, right_novel$summary$contactRight)
  expect_equal(right_novel$summary$contactFamiliar, right_novel$summary$contactLeft)

  left_novel <- compute_socp_metrics(tracking, "R", fps = fps)
  expect_equal(left_novel$summary$contactNovel, left_novel$summary$contactLeft)

  missing <- compute_socp_metrics(tracking, NA_character_, fps = fps)
  expect_true(is.na(missing$summary$contactNovel))
  expect_equal(missing$summary$contactLeft, 2)
})

test_that("no contact yields censored NA latency", {
  n <- 20
  tracking <- socp_tracking(
    n = n, fps = 10, nose = c(0, 0), bodycentre = c(0, -5),
    socl = c(-200, 0), socr = c(200, 0)
  )
  result <- compute_socp_metrics(tracking, "R", fps = 10)
  expect_equal(result$summary$contactLeft, 0)
  expect_true(is.na(result$summary$latencyLeft))
  expect_false(is.infinite(result$summary$latencyLeft))
  expect_equal(result$summary$frequencyLeft, 0L)
})

test_that("bout structure matches the constructed trajectory", {
  n <- 20
  fps <- 10
  near <- c(rep(TRUE, 4), rep(FALSE, 6), rep(TRUE, 5), rep(FALSE, 5))
  nose_x <- ifelse(near, 2, 100)
  tracking <- socp_tracking(
    n = n, fps = fps,
    nose = cbind(x = nose_x, y = rep(0, n)),
    bodycentre = cbind(x = rep(-5, n), y = rep(0, n)),
    socl = c(2, 0), socr = c(300, 0)
  )
  result <- compute_socp_metrics(tracking, "R", fps = fps)
  expect_equal(result$summary$frequencyLeft, 2L)
  expect_equal(result$summary$contactLeft, 0.9)
  expect_equal(result$summary$latencyLeft, 0)
  expect_equal(result$summary$entriesLeft, 1L)
})

test_that("the body exclusion distance rejects sitting on the stimulus", {
  n <- 20
  tracking <- socp_tracking(
    n = n, fps = 10, nose = c(2, 0), bodycentre = c(2, 0.5),
    socl = c(2, 0), socr = c(300, 0)
  )
  expect_equal(compute_socp_metrics(tracking, "R", fps = 10)$summary$contactLeft, 0)
})

test_that("optional orientation gating is off by default", {
  n <- 20
  # bodycentre -> nose points straight at the stimulus, which the NOR angle
  # convention rejects. SocP keeps its established definition by default.
  tracking <- socp_tracking(
    n = n, fps = 10, nose = c(2, 0), bodycentre = c(2, 4),
    socl = c(2, 0), socr = c(300, 0)
  )
  expect_gt(compute_socp_metrics(tracking, "R", fps = 10)$summary$contactLeft, 0)
  expect_equal(
    compute_socp_metrics(tracking, "R", fps = 10, require_orientation = TRUE)$summary$contactLeft,
    0
  )
})
