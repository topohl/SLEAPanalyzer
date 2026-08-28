# Zone membership, entries and valid time on synthetic trajectories.

zone_tracking <- function(x, y, fps = 10, zones = NULL) {
  n <- length(x)
  points <- list(bodycentre = cbind(x = x, y = y))
  tracking <- make_tracking(frames = 0:(n - 1), fps = fps, points = points)
  if (is.null(zones)) {
    zones <- list(
      center = data.frame(x = c(3, 7, 7, 3), y = c(3, 3, 7, 7)),
      periphery = data.frame(x = c(3, 7, 7, 3), y = c(3, 3, 7, 7))
    )
  }
  tracking$zones <- zones
  tracking
}

test_that("unobserved frames belong to no zone and to no inverted zone", {
  x <- c(5, NA, 5, 0, 0)
  y <- c(5, NA, 5, 0, 0)
  tracking <- zone_tracking(x, y)

  inside <- IsInZone(tracking, "bodycentre", "center")
  expect_equal(inside, c(TRUE, FALSE, TRUE, FALSE, FALSE))

  # The decisive case: sp reports an untracked frame as outside every polygon,
  # so inverting before masking credited every dropout to the complementary
  # zone. Frame 2 must belong to neither.
  outside <- IsInZone(tracking, "bodycentre", "center", invert = TRUE)
  expect_equal(outside, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_false(any(inside & outside))
  expect_false(inside[2] || outside[2])
})

test_that("ZoneReport separates entries from legacy transitions", {
  # Two visits to the centre, starting and ending outside.
  x <- c(0, 0, 5, 5, 0, 0, 5, 5, 0, 0)
  tracking <- zone_tracking(x, rep(5, length(x)))
  tracking <- CalculateMovement(tracking, movement_cutoff = 1, integration_period = 0)
  tracking$integration_period <- 0

  report <- ZoneReport(tracking, "bodycentre", "center", zone.name = "center")

  expect_equal(report$center.entries, 2L)
  # Legacy definition: onsets plus offsets.
  expect_equal(report$center.transitions, 4L)
  expect_equal(report$center.valid.time, 1.0)
})

test_that("a visit in progress at the first frame is not counted as an entry", {
  x <- c(5, 5, 0, 0, 5, 5)
  tracking <- zone_tracking(x, rep(5, length(x)))
  tracking <- CalculateMovement(tracking, movement_cutoff = 1, integration_period = 0)
  tracking$integration_period <- 0

  report <- ZoneReport(tracking, "bodycentre", "center", zone.name = "center")
  # One observed onset (frame 5); the opening visit has no observed onset.
  expect_equal(report$center.entries, 1L)
  expect_equal(report$center.transitions, 2L)
})

test_that("zone occupancy time excludes unobserved frames", {
  x <- c(5, 5, NA, NA, 5, 5, 0, 0, 0, 0)
  tracking <- zone_tracking(x, rep(5, length(x)))
  tracking <- CalculateMovement(tracking, movement_cutoff = 1, integration_period = 0)
  tracking$integration_period <- 0

  report <- ZoneReport(tracking, "bodycentre", "center", zone.name = "center")
  # Four observed frames in the centre out of eight observed frames overall.
  expect_equal(report$center.total.time, 0.4)
  expect_equal(report$center.valid.time, 0.8)
})

test_that("a dropout inside a visit does not fabricate a zone entry", {
  x <- c(0, 5, 5, NA, 5, 5, 0)
  tracking <- zone_tracking(x, rep(5, length(x)))
  tracking <- CalculateMovement(tracking, movement_cutoff = 1, integration_period = 0)
  tracking$integration_period <- 0

  report <- ZoneReport(tracking, "bodycentre", "center", zone.name = "center")
  expect_equal(report$center.entries, 1L)
})
