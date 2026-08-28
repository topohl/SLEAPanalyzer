# Synthetic elevated-plus-maze runs with analytically known results.
#
# Maze layout in calibrated units, a plus centred on (50, 50) with arms 20
# wide and reaching to 0 and 100:
#
#             open.top      y 80..100
#   closed.left  center  closed.right    y 40..60
#             open.bottom   y 0..20

epm_zones <- function() {
  list(
    arena = data.frame(
      x = c(40, 60, 60, 100, 100, 60, 60, 40, 40, 0, 0, 40),
      y = c(100, 100, 60, 60, 40, 40, 0, 0, 40, 40, 60, 60)
    ),
    center = data.frame(x = c(40, 60, 60, 40), y = c(40, 40, 60, 60)),
    open.top = data.frame(x = c(40, 60, 60, 40), y = c(60, 60, 100, 100)),
    open.bottom = data.frame(x = c(40, 60, 60, 40), y = c(0, 0, 40, 40)),
    closed.left = data.frame(x = c(0, 40, 40, 0), y = c(40, 40, 60, 60)),
    closed.right = data.frame(x = c(60, 100, 100, 60), y = c(40, 40, 60, 60))
  )
}

epm_tracking <- function(body_x, body_y, fps = 10,
                         head_x = NULL, head_y = NULL,
                         neck_x = NULL, neck_y = NULL) {
  n <- length(body_x)
  if (is.null(head_x)) head_x <- body_x
  if (is.null(head_y)) head_y <- body_y
  if (is.null(neck_x)) neck_x <- body_x
  if (is.null(neck_y)) neck_y <- body_y
  points <- list(
    bodycentre = cbind(x = body_x, y = body_y),
    headcentre = cbind(x = head_x, y = head_y),
    neck = cbind(x = neck_x, y = neck_y)
  )
  tracking <- make_tracking(frames = 0:(n - 1), fps = fps, points = points,
                            distance_unit = "cm")
  tracking$zones <- epm_zones()
  tracking
}

test_that("EPMAnalysis runs on synthetic tracking and assigns arms correctly", {
  # 4 frames in the closed left arm, 4 in the centre, 4 in the open top arm.
  body_x <- c(rep(20, 4), rep(50, 4), rep(50, 4))
  body_y <- c(rep(50, 4), rep(50, 4), rep(80, 4))
  tracking <- epm_tracking(body_x, body_y, fps = 4)

  result <- EPMAnalysis(
    tracking, movement_cutoff = 5, integration_period = 0,
    points = "bodycentre", nosedips = TRUE
  )

  expect_equal(result$Report$bodycentre.closed.total.time, 1.0)
  expect_equal(result$Report$bodycentre.center.total.time, 1.0)
  expect_equal(result$Report$bodycentre.open.total.time, 1.0)
  expect_equal(result$Report$bodycentre.closed.left.total.time, 1.0)
  expect_equal(result$Report$bodycentre.closed.right.total.time, 0)
  expect_equal(result$Report$bodycentre.open.top.total.time, 1.0)
  expect_equal(result$Report$bodycentre.open.bottom.total.time, 0)
  expect_equal(result$Report$bodycentre.total.time, 3.0)
  expect_equal(result$Report$bodycentre.valid.time, 3.0)
})

test_that("EPM entries count observed onsets and transitions stay available", {
  # Out, into the open arm, out, into the open arm again.
  body_x <- rep(50, 12)
  body_y <- c(rep(50, 3), rep(80, 3), rep(50, 3), rep(80, 3))
  tracking <- epm_tracking(body_x, body_y, fps = 3)

  result <- EPMAnalysis(
    tracking, movement_cutoff = 5, integration_period = 0, points = "bodycentre"
  )

  expect_equal(result$Report$bodycentre.open.entries, 2L)
  # Legacy formula counts onsets plus offsets.
  expect_equal(result$Report$bodycentre.open.transitions, 3L)
})

test_that("EPM excludes unobserved frames from arm times and valid time", {
  body_x <- c(rep(20, 4), rep(NA_real_, 4), rep(50, 4))
  body_y <- c(rep(50, 4), rep(NA_real_, 4), rep(80, 4))
  tracking <- epm_tracking(body_x, body_y, fps = 4)

  result <- EPMAnalysis(
    tracking, movement_cutoff = 5, integration_period = 0, points = "bodycentre"
  )

  expect_equal(result$Report$bodycentre.total.time, 3.0)
  expect_equal(result$Report$bodycentre.valid.time, 2.0)
  expect_equal(result$Report$bodycentre.valid.fraction, 2 / 3)
  expect_equal(result$Report$bodycentre.closed.total.time, 1.0)
  expect_equal(result$Report$bodycentre.open.total.time, 1.0)
  # The dropout must not be credited to the centre, which is what the
  # complementary-zone bug would have done.
  expect_equal(result$Report$bodycentre.center.total.time, 0)
})

test_that("nose dips are counted as whole onsets, never half events", {
  n <- 12
  # The body stays on the open arm; the head leaves the maze outline on
  # frames 4-6 and again on frames 10-12, the second dip still in progress
  # at the last frame.
  body_x <- rep(50, n)
  body_y <- rep(80, n)
  head_x <- rep(50, n)
  head_y <- rep(80, n)
  head_x[4:6] <- 150   # off the maze entirely
  head_x[10:12] <- 150
  tracking <- epm_tracking(body_x, body_y, fps = 3,
                           head_x = head_x, head_y = head_y)

  result <- EPMAnalysis(
    tracking, movement_cutoff = 5, integration_period = 0,
    points = "bodycentre", nosedips = TRUE
  )

  expect_equal(result$Report$nose.dip, 2L)
  # The legacy CalculateTransitions(...) / 2 returned 1.5 here because the
  # trailing dip contributed an onset without a matching offset.
  expect_true(result$Report$nose.dip == round(result$Report$nose.dip))
  expect_equal(result$Report$nose.dip.valid.time, 4.0)
})

test_that("a nose dip requires the neck outside the closed arms", {
  n <- 9
  # Head off the maze but neck in the closed left arm: not a nose dip.
  body_x <- rep(20, n)
  body_y <- rep(50, n)
  head_x <- rep(-50, n)
  head_y <- rep(50, n)
  tracking <- epm_tracking(body_x, body_y, fps = 3,
                           head_x = head_x, head_y = head_y,
                           neck_x = rep(20, n), neck_y = rep(50, n))

  result <- EPMAnalysis(
    tracking, movement_cutoff = 5, integration_period = 0,
    points = "bodycentre", nosedips = TRUE
  )
  expect_equal(result$Report$nose.dip, 0L)
})

test_that("an unobserved head does not register as a nose dip", {
  n <- 9
  body_x <- rep(50, n)
  body_y <- rep(80, n)
  head_x <- rep(50, n)
  head_y <- rep(80, n)
  head_x[4:6] <- NA_real_
  head_y[4:6] <- NA_real_
  tracking <- epm_tracking(body_x, body_y, fps = 3,
                           head_x = head_x, head_y = head_y)

  result <- EPMAnalysis(
    tracking, movement_cutoff = 5, integration_period = 0,
    points = "bodycentre", nosedips = TRUE
  )
  # IsInZone() reports FALSE for a dropout, and the dip test negates the head
  # term, so an unmasked implementation would score the dropout as a dip.
  expect_equal(result$Report$nose.dip, 0L)
  expect_equal(result$Report$nose.dip.valid.time, 2.0)
})

test_that("EPMAnalysis warns when nose-dip landmarks or zones are missing", {
  tracking <- epm_tracking(rep(50, 6), rep(50, 6), fps = 3)
  tracking$data$neck <- NULL
  expect_warning(
    EPMAnalysis(tracking, movement_cutoff = 5, integration_period = 0,
                points = "bodycentre", nosedips = TRUE),
    "nosedip analysis"
  )
})

test_that("EPM distance and speed are reported in calibrated units", {
  # Ten frames moving 2 cm per frame at 10 fps is 20 cm/s and 18 cm total
  # (nine real displacement intervals).
  body_x <- seq(40, by = 2, length.out = 10)
  tracking <- epm_tracking(body_x, rep(50, 10), fps = 10)

  result <- EPMAnalysis(
    tracking, movement_cutoff = 5, integration_period = 0, points = "bodycentre"
  )
  expect_equal(result$Report$bodycentre.raw.distance, 18)
  expect_equal(result$Report$bodycentre.raw.speed, 20)
})
