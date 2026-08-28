test_that("CleanTrackingData leaves long gaps missing instead of fabricating them", {
  n <- 40
  fps <- 10
  points <- list(bodycentre = cbind(x = as.numeric(1:n), y = rep(0, n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = fps, points = points)
  # A one-second dropout at 10 fps.
  tracking$data$bodycentre$likelihood[15:24] <- 0.1

  cleaned <- CleanTrackingData(tracking, likelihoodcutoff = 0.9, max.gap.s = 0.3)
  expect_true(all(is.na(cleaned$data$bodycentre$x[15:24])))
  expect_equal(as.character(cleaned$data$bodycentre$status[15:24]), rep("invalid", 10))
  expect_equal(as.character(cleaned$data$bodycentre$status[1]), "observed")
})

test_that("CleanTrackingData interpolates short gaps and labels them", {
  n <- 20
  points <- list(bodycentre = cbind(x = as.numeric(1:n), y = rep(0, n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  tracking$data$bodycentre$likelihood[10:11] <- 0.1

  cleaned <- CleanTrackingData(tracking, likelihoodcutoff = 0.9, max.gap.s = 0.3)
  expect_equal(cleaned$data$bodycentre$x[10:11], c(10, 11))
  expect_equal(as.character(cleaned$data$bodycentre$status[10:11]),
               c("interpolated", "interpolated"))
})

test_that("CleanTrackingData never forward-fills leading or trailing gaps", {
  n <- 20
  points <- list(bodycentre = cbind(x = as.numeric(1:n), y = rep(0, n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  tracking$data$bodycentre$likelihood[1:3] <- 0.1
  tracking$data$bodycentre$likelihood[18:20] <- 0.1

  cleaned <- CleanTrackingData(tracking, likelihoodcutoff = 0.9, max.gap.s = 10)
  expect_true(all(is.na(cleaned$data$bodycentre$x[1:3])))
  expect_true(all(is.na(cleaned$data$bodycentre$x[18:20])))
})

test_that("CleanTrackingData preserves the likelihood column", {
  n <- 20
  points <- list(bodycentre = cbind(x = as.numeric(1:n), y = rep(0, n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  tracking$data$bodycentre$likelihood[10] <- 0.1

  cleaned <- CleanTrackingData(tracking, likelihoodcutoff = 0.9, max.gap.s = 0.3)
  # The old implementation interpolated the whole data frame, overwriting the
  # confidence value that recorded why the frame had been rejected.
  expect_equal(cleaned$data$bodycentre$likelihood[10], 0.1)
})

test_that("CleanTrackingData warns when asked for unbounded interpolation", {
  n <- 10
  points <- list(bodycentre = cbind(x = as.numeric(1:n), y = rep(0, n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  expect_warning(
    CleanTrackingData(tracking, likelihoodcutoff = 0.9, max.gap.s = Inf),
    "unlimited length"
  )
  expect_error(
    CleanTrackingData(tracking, likelihoodcutoff = 0.9, max.gap.s = -1),
    "max.gap.s"
  )
})

test_that("CleanTrackingData rejects implausible jumps without interpolating across them", {
  n <- 20
  x <- as.numeric(1:n)
  x[10] <- 500
  points <- list(bodycentre = cbind(x = x, y = rep(0, n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)

  cleaned <- CleanTrackingData(
    tracking, likelihoodcutoff = 0.5, maxdelta = 5, max.gap.s = 0.3
  )
  # The teleport is rejected and replaced by the interpolated value.
  expect_equal(cleaned$data$bodycentre$x[10], 10)
  expect_equal(as.character(cleaned$data$bodycentre$status[10]), "interpolated")
})

test_that("zone membership includes points exactly on the boundary", {
  n <- 5
  zone <- data.frame(x = c(0, 10, 10, 0), y = c(0, 0, 10, 10))
  points <- list(
    # Frame 1 is on an edge, frame 2 on a vertex, frame 3 strictly inside,
    # frames 4 and 5 outside.
    bodycentre = cbind(x = c(5, 0, 5, -1, 20), y = c(0, 0, 5, 5, 20))
  )
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  tracking$zones <- list(arena = zone)

  # The legacy test for == 1 reported FALSE for the edge and vertex frames,
  # so an animal on a shared zone boundary belonged to no zone at all.
  expect_equal(
    IsInZone(tracking, "bodycentre", "arena"),
    c(TRUE, TRUE, TRUE, FALSE, FALSE)
  )
})
