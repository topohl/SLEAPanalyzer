testthat::test_that("NOR uses the corresponding object angle and frame-derived timing", {
  n <- 60
  points <- square_points(n)
  points$bodycentre <- cbind(x = rep(-1, n), y = rep(0, n))
  points$nose <- cbind(x = rep(0, n), y = rep(0, n))
  points$objL <- cbind(x = rep(2, n), y = rep(0, n))
  points$objR <- cbind(x = rep(-3, n), y = rep(0, n))
  result <- compute_nor_metrics(
    make_tracking(frames = 0:(n - 1), points = points, distance_unit = "cm"),
    "R", fps = 30
  )

  testthat::expect_equal(result$summary$contactLeft, 2)
  testthat::expect_equal(result$summary$contactRight, 0)
  testthat::expect_equal(result$summary$contactNov, result$summary$contactLeft)
  testthat::expect_equal(result$summary$frequencyL, 1)
  testthat::expect_equal(result$summary$latencyLeft, 0)
  testthat::expect_equal(result$summary$totalTime, 2)
  testthat::expect_equal(result$angles$right, rep(0, n))
  # Every pre-v2 column is still emitted; v2 only appends columns.
  legacy_columns <- c(
    "contactLeft", "contactRight", "contactNov", "contactFam", "proxLeft",
    "proxRight", "proxNov", "proxFam", "proxLeftAngle", "proxRightAngle",
    "latency", "latencyLeft", "latencyRight", "frequencyL", "frequencyR",
    "totalTime", "novelLoc"
  )
  testthat::expect_true(all(legacy_columns %in% names(result$summary)))
  testthat::expect_named(
    result$summary,
    c(
      "contactLeft", "contactRight", "contactNov", "contactFam", "proxLeft",
      "proxRight", "proxNov", "proxFam", "proxLeftAngle", "proxRightAngle",
      "latency", "latencyLeft", "latencyRight", "frequencyL", "frequencyR",
      "meanBoutLeft", "meanBoutRight", "entriesLeft", "entriesRight",
      "totalTime", "validTime", "validFraction", "contactGeometry", "novelLoc"
    )
  )
})

testthat::test_that("NOR no-contact latency is NA and separated visits count as entries", {
  n <- 10
  points <- square_points(n)
  points$bodycentre <- cbind(x = rep(-10, n), y = 0)
  points$nose <- cbind(x = c(0, 0, 20, 20, 0, 0, 20, 20, 20, 20), y = 0)
  points$objL <- cbind(x = rep(2, n), y = 0)
  points$objR <- cbind(x = rep(100, n), y = 0)
  result <- compute_nor_metrics(
    make_tracking(points = points, fps = 10, distance_unit = "cm"), "R", fps = 10
  )
  testthat::expect_equal(result$summary$frequencyL, 2)
  testthat::expect_true(is.na(result$summary$latencyRight))

  points$nose[,] <- 50
  absent <- compute_nor_metrics(
    make_tracking(points = points, fps = 10, distance_unit = "cm"), "R", fps = 10
  )
  testthat::expect_true(is.na(absent$summary$latency))
  testthat::expect_false(is.infinite(absent$summary$latency))
})

testthat::test_that("SocP contact, remapping, entries, and missing metadata are safe", {
  n <- 8
  points <- square_points(n)
  points$bodycentre <- cbind(x = rep(-5, n), y = 0)
  points$nose <- cbind(x = c(0, 0, 20, 20, 0, 0, 20, 20), y = 0)
  points$socl <- cbind(x = rep(2, n), y = 0)
  points$socr <- cbind(x = rep(100, n), y = 0)
  tracking <- make_tracking(frames = 0:(n - 1), points = points)

  result <- compute_socp_metrics(tracking, "R", fps = 2)
  testthat::expect_equal(result$summary$contactLeft, 2)
  testthat::expect_equal(result$summary$contactNovel, result$summary$contactLeft)
  testthat::expect_equal(result$summary$frequencyLeft, 2)
  testthat::expect_true(is.na(result$summary$latencyRight))

  missing <- compute_socp_metrics(tracking, NA_character_, fps = 2)
  testthat::expect_true(is.na(missing$summary$contactNovel))
  testthat::expect_equal(missing$summary$contactLeft, 2)
})
