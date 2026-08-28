# Synthetic NOR trajectories with analytically known results.
#
# The arena landmarks come from square_points(); only the animal and object
# positions carry the behavioral signal.

nor_tracking <- function(n = 60, fps = 30, nose, bodycentre, objL, objR,
                         unit = "cm") {
  points <- square_points(n)
  rep_xy <- function(xy) {
    if (is.matrix(xy)) return(xy)
    cbind(x = rep(xy[1], n), y = rep(xy[2], n))
  }
  points$nose <- rep_xy(nose)
  points$bodycentre <- rep_xy(bodycentre)
  points$objL <- rep_xy(objL)
  points$objR <- rep_xy(objR)
  make_tracking(frames = 0:(n - 1), fps = fps, points = points, distance_unit = unit)
}

test_that("the same contact detector is applied to both objects", {
  # The nose sits exactly halfway between two objects, 4.2 cm from each, so
  # the true behavior toward the two objects is identical by construction.
  # 4.2 is inside the 9x7 box (half-width 4.5) but outside the 4 cm radius.
  n <- 30
  tracking <- nor_tracking(
    n = n, fps = 30,
    nose = c(0, 0), bodycentre = c(0, -3),
    objL = c(-4.2, 0), objR = c(4.2, 0)
  )

  radial <- compute_nor_metrics(tracking, "R", fps = 30, contact_geometry = "radial")
  expect_equal(radial$summary$contactLeft, radial$summary$contactRight)
  expect_equal(radial$summary$contactNov, radial$summary$contactFam)
  expect_equal(radial$summary$contactNov, 0)

  box <- compute_nor_metrics(tracking, "R", fps = 30, contact_geometry = "box")
  expect_equal(box$summary$contactLeft, box$summary$contactRight)
  expect_equal(box$summary$contactNov, box$summary$contactFam)
  expect_equal(box$summary$contactNov, 1)
})

test_that("the legacy asymmetric detector manufactures a discrimination index", {
  n <- 30
  tracking <- nor_tracking(
    n = n, fps = 30,
    nose = c(0, 0), bodycentre = c(0, -3),
    objL = c(-4.2, 0), objR = c(4.2, 0)
  )
  legacy <- expect_warning(
    compute_nor_metrics(tracking, "R", fps = 30, contact_geometry = "legacy_asymmetric"),
    "biased"
  )
  # Identical behavior toward both objects, yet the novel object scores a full
  # second of contact and the familiar object scores none. This is the defect
  # the symmetric detector removes.
  expect_equal(legacy$summary$contactNov, 1)
  expect_equal(legacy$summary$contactFam, 0)
})

test_that("novel mapping follows the metadata for a left novel object", {
  n <- 30
  # Nose parked on the left object only.
  tracking <- nor_tracking(
    n = n, fps = 30,
    nose = c(-2, 0), bodycentre = c(-2, -3),
    objL = c(-2, 2), objR = c(40, 0)
  )
  # Metadata "R" maps the left side to novel in this repository's convention.
  result <- compute_nor_metrics(tracking, "R", fps = 30)
  expect_equal(result$summary$novelLoc, "R")
  expect_gt(result$summary$contactLeft, 0)
  expect_equal(result$summary$contactRight, 0)
  expect_equal(result$summary$contactNov, result$summary$contactLeft)
  expect_equal(result$summary$contactFam, result$summary$contactRight)
})

test_that("novel mapping follows the metadata for a right novel object", {
  n <- 30
  tracking <- nor_tracking(
    n = n, fps = 30,
    nose = c(-2, 0), bodycentre = c(-2, -3),
    objL = c(-2, 2), objR = c(40, 0)
  )
  result <- compute_nor_metrics(tracking, "L", fps = 30)
  expect_equal(result$summary$contactNov, result$summary$contactRight)
  expect_equal(result$summary$contactFam, result$summary$contactLeft)
})

test_that("contact is measured even when novel-location metadata is missing", {
  n <- 30
  tracking <- nor_tracking(
    n = n, fps = 30,
    nose = c(-2, 0), bodycentre = c(-2, -3),
    objL = c(-2, 2), objR = c(40, 0)
  )
  result <- compute_nor_metrics(tracking, NA_character_, fps = 30)
  # Side-specific measurements survive; only the novel/familiar mapping is NA.
  expect_gt(result$summary$contactLeft, 0)
  expect_true(is.na(result$summary$contactNov))
  expect_true(is.na(result$summary$contactFam))
  expect_true(is.na(result$summary$novelLoc))
})

test_that("no contact yields censored NA latency, never Inf", {
  n <- 30
  tracking <- nor_tracking(
    n = n, fps = 30,
    nose = c(0, 0), bodycentre = c(0, -3),
    objL = c(-40, 0), objR = c(40, 0)
  )
  result <- compute_nor_metrics(tracking, "R", fps = 30)
  expect_equal(result$summary$contactLeft, 0)
  expect_true(is.na(result$summary$latency))
  expect_true(is.na(result$summary$latencyLeft))
  expect_false(is.infinite(result$summary$latency))
  expect_equal(result$summary$frequencyL, 0L)
})

test_that("bout structure and latency match the constructed trajectory", {
  n <- 20
  fps <- 10
  # Frames 1-4 and 11-15 (1-based) are on the object; the rest are far away.
  on_object <- c(rep(FALSE, 0), rep(TRUE, 4), rep(FALSE, 6), rep(TRUE, 5), rep(FALSE, 5))
  nose_x <- ifelse(on_object, -2, 40)
  tracking <- nor_tracking(
    n = n, fps = fps,
    nose = cbind(x = nose_x, y = rep(0, n)),
    bodycentre = cbind(x = nose_x, y = rep(-3, n)),
    objL = c(-2, 2), objR = c(200, 0)
  )
  result <- compute_nor_metrics(tracking, "R", fps = fps)
  expect_equal(result$summary$frequencyL, 2L)
  expect_equal(result$summary$contactLeft, 9 / fps)
  expect_equal(result$summary$latencyLeft, 0)
  expect_equal(result$summary$entriesLeft, 1L)
  expect_equal(result$summary$meanBoutLeft, 4.5 / fps)
})

test_that("minimum bout duration is applied in seconds", {
  n <- 20
  fps <- 10
  on_object <- c(rep(TRUE, 2), rep(FALSE, 8), rep(TRUE, 5), rep(FALSE, 5))
  nose_x <- ifelse(on_object, -2, 40)
  tracking <- nor_tracking(
    n = n, fps = fps,
    nose = cbind(x = nose_x, y = rep(0, n)),
    bodycentre = cbind(x = nose_x, y = rep(-3, n)),
    objL = c(-2, 2), objR = c(200, 0)
  )
  # A 0.3 s minimum drops the two-frame (0.2 s) bout but keeps the five-frame one.
  result <- compute_nor_metrics(tracking, "R", fps = fps, min_bout_s = 0.3)
  expect_equal(result$summary$frequencyL, 1L)
  expect_equal(result$summary$contactLeft, 0.5)
  expect_equal(result$summary$latencyLeft, 1.0)
})

test_that("missing frames reduce valid time and are not scored as contact", {
  n <- 30
  fps <- 10
  tracking <- nor_tracking(
    n = n, fps = fps,
    nose = c(-2, 0), bodycentre = c(-2, -3),
    objL = c(-2, 2), objR = c(200, 0)
  )
  # Ten frames of the nose were never tracked.
  tracking$data$nose$x[11:20] <- NA_real_
  tracking$data$nose$y[11:20] <- NA_real_

  result <- compute_nor_metrics(tracking, "R", fps = fps)
  expect_equal(result$summary$totalTime, 3)
  expect_equal(result$summary$validTime, 2)
  expect_equal(result$summary$validFraction, 2 / 3)
  # Contact is credited only for the 20 observed frames, not all 30.
  expect_equal(result$summary$contactLeft, 2)
  # The dropout splits the visit into two observed bouts.
  expect_equal(result$summary$frequencyL, 2L)
})

test_that("a bridged dropout rejoins one visit without inflating its duration", {
  n <- 30
  fps <- 10
  tracking <- nor_tracking(
    n = n, fps = fps,
    nose = c(-2, 0), bodycentre = c(-2, -3),
    objL = c(-2, 2), objR = c(200, 0)
  )
  tracking$data$nose$x[15:16] <- NA_real_
  tracking$data$nose$y[15:16] <- NA_real_

  result <- compute_nor_metrics(tracking, "R", fps = fps, max_gap_s = 0.5)
  expect_equal(result$summary$frequencyL, 1L)
  # 28 observed frames, not 30: bridging joins bouts but never invents data.
  expect_equal(result$summary$contactLeft, 2.8)
  expect_equal(result$summary$validTime, 2.8)
})

test_that("orientation away from the object suppresses contact", {
  n <- 30
  # bodycentre -> nose points straight at the object, so the legacy angle
  # convention gives 0 degrees and the >= 70 criterion rejects the frame.
  facing_away <- nor_tracking(
    n = n, fps = 30,
    nose = c(-2, 0), bodycentre = c(-2, 4), objL = c(-2, 2), objR = c(200, 0)
  )
  expect_equal(compute_nor_metrics(facing_away, "R", fps = 30)$summary$contactLeft, 0)

  facing_object <- nor_tracking(
    n = n, fps = 30,
    nose = c(-2, 0), bodycentre = c(-2, -3), objL = c(-2, 2), objR = c(200, 0)
  )
  expect_gt(compute_nor_metrics(facing_object, "R", fps = 30)$summary$contactLeft, 0)
})

test_that("the body exclusion distance rejects the animal sitting on the object", {
  n <- 30
  # bodycentre within 1 cm of the object: the nose position is not evidence of
  # investigation because the whole animal is on top of the object.
  tracking <- nor_tracking(
    n = n, fps = 30,
    nose = c(-2, 0), bodycentre = c(-2, 1.5), objL = c(-2, 2), objR = c(200, 0)
  )
  expect_equal(compute_nor_metrics(tracking, "R", fps = 30)$summary$contactLeft, 0)
})

test_that("NOR rejects pixel coordinates when thresholds are in cm", {
  tracking <- nor_tracking(
    n = 10, fps = 10,
    nose = c(0, 0), bodycentre = c(0, -3), objL = c(-2, 0), objR = c(40, 0),
    unit = "pixel"
  )
  expect_error(compute_nor_metrics(tracking, "R", fps = 10), "does not match")
})

test_that("NOR rejects an fps that disagrees with the tracking data", {
  tracking <- nor_tracking(
    n = 10, fps = 10,
    nose = c(0, 0), bodycentre = c(0, -3), objL = c(-2, 0), objR = c(40, 0)
  )
  expect_error(compute_nor_metrics(tracking, "R", fps = 30), "fps does not match")
})
