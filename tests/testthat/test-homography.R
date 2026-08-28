# Arena rectification against synthetic cameras with known ground truth.

# Corner order matches the repository convention: tl, tr, br, bl.
unit_square_corners <- function() {
  data.frame(x = c(0, 1, 1, 0), y = c(1, 1, 0, 0))
}

#' Project world points through a known homography to simulate a camera.
project_with <- function(H, x, y) {
  w <- H[3, 1] * x + H[3, 2] * y + H[3, 3]
  list(
    x = (H[1, 1] * x + H[1, 2] * y + H[1, 3]) / w,
    y = (H[2, 1] * x + H[2, 2] * y + H[2, 3]) / w
  )
}

test_that("an identity mapping is recovered exactly", {
  corners <- rectangular_arena_target(50, 50)
  calibration <- arena_calibration(corners, 50, 50)
  expect_lt(calibration$max_reprojection_error, 1e-9)

  mapped <- rectify_coordinates(calibration, c(0, 25, 50), c(0, 25, 50))
  expect_equal(mapped$x, c(0, 25, 50))
  expect_equal(mapped$y, c(0, 25, 50))
})

test_that("translation is removed and the origin lands on the arena corner", {
  # Camera sees the arena offset by (1000, 500) pixels.
  corners <- data.frame(x = c(1000, 1500, 1500, 1000), y = c(1000, 1000, 500, 500))
  calibration <- arena_calibration(corners, 50, 50)

  bottom_left <- rectify_coordinates(calibration, 1000, 500)
  expect_equal(bottom_left$x, 0)
  expect_equal(bottom_left$y, 0)

  top_right <- rectify_coordinates(calibration, 1500, 1000)
  expect_equal(top_right$x, 50)
  expect_equal(top_right$y, 50)
})

test_that("anisotropic scaling is corrected independently per axis", {
  # 10 px/cm horizontally but 4 px/cm vertically: a single scalar scale
  # cannot represent this, which is why scalar calibration is insufficient.
  corners <- data.frame(x = c(0, 500, 500, 0), y = c(200, 200, 0, 0))
  calibration <- arena_calibration(corners, 50, 50)

  centre <- rectify_coordinates(calibration, 250, 100)
  expect_equal(centre$x, 25)
  expect_equal(centre$y, 25)

  # A 100 px horizontal step is 10 cm; a 100 px vertical step is 25 cm.
  a <- rectify_coordinates(calibration, 100, 100)
  b <- rectify_coordinates(calibration, 200, 100)
  expect_equal(b$x - a$x, 10)
  c1 <- rectify_coordinates(calibration, 100, 50)
  c2 <- rectify_coordinates(calibration, 100, 150)
  expect_equal(c2$y - c1$y, 25)
})

test_that("rotation is removed and distances are preserved", {
  theta <- pi / 6
  world <- rectangular_arena_target(40, 40)
  rotated <- data.frame(
    x = world$x * cos(theta) - world$y * sin(theta) + 300,
    y = world$x * sin(theta) + world$y * cos(theta) + 120
  )
  calibration <- arena_calibration(rotated, 40, 40)
  expect_lt(calibration$max_reprojection_error, 1e-9)

  # Two world points 10 cm apart must remain 10 cm apart after rectification.
  p <- c(10, 20)
  q <- c(20, 20)
  to_image <- function(pt) c(
    pt[1] * cos(theta) - pt[2] * sin(theta) + 300,
    pt[1] * sin(theta) + pt[2] * cos(theta) + 120
  )
  ip <- to_image(p); iq <- to_image(q)
  rp <- rectify_coordinates(calibration, ip[1], ip[2])
  rq <- rectify_coordinates(calibration, iq[1], iq[2])
  expect_equal(sqrt((rp$x - rq$x)^2 + (rp$y - rq$y)^2), 10)
})

test_that("perspective distortion is removed and recovers true distances", {
  # A genuine projective transform: the last row is non-zero, so the
  # pixels-per-centimetre ratio varies across the image.
  H <- matrix(c(
    8, 0.6, 100,
    0.4, 7.5, 80,
    0.002, 0.004, 1
  ), nrow = 3, byrow = TRUE)

  world <- rectangular_arena_target(60, 40)
  image_corners <- project_with(H, world$x, world$y)
  calibration <- arena_calibration(
    data.frame(x = image_corners$x, y = image_corners$y), 60, 40
  )
  expect_lt(calibration$max_reprojection_error, 1e-8)
  expect_gt(calibration$perspective_index, 0)

  # Sample a grid of world points, project them, rectify them back and check
  # the round trip. A scalar scale cannot do this.
  grid <- expand.grid(x = seq(2, 58, by = 8), y = seq(2, 38, by = 6))
  imaged <- project_with(H, grid$x, grid$y)
  recovered <- rectify_coordinates(calibration, imaged$x, imaged$y)
  expect_equal(recovered$x, grid$x, tolerance = 1e-8)
  expect_equal(recovered$y, grid$y, tolerance = 1e-8)

  # Ground-truth distance between two known world points.
  p <- c(10, 10); q <- c(40, 30)
  true_distance <- sqrt(sum((p - q)^2))
  ip <- project_with(H, p[1], p[2])
  iq <- project_with(H, q[1], q[2])
  rp <- rectify_coordinates(calibration, ip$x, ip$y)
  rq <- rectify_coordinates(calibration, iq$x, iq$y)
  expect_equal(sqrt((rp$x - rq$x)^2 + (rp$y - rq$y)^2), true_distance, tolerance = 1e-8)
})

test_that("a scalar scale is measurably wrong under perspective", {
  # Demonstrates why the homography matters rather than merely asserting it.
  H <- matrix(c(
    8, 0, 100,
    0, 8, 80,
    0.004, 0.004, 1
  ), nrow = 3, byrow = TRUE)
  world <- rectangular_arena_target(60, 60)
  image_corners <- project_with(H, world$x, world$y)
  calibration <- arena_calibration(
    data.frame(x = image_corners$x, y = image_corners$y), 60, 60
  )

  # A 10 cm rod near the near edge and the same rod near the far edge.
  rod <- function(x0, y0) {
    a <- project_with(H, x0, y0)
    b <- project_with(H, x0 + 10, y0)
    c(pixels = sqrt((a$x - b$x)^2 + (a$y - b$y)^2))
  }
  near_px <- rod(5, 5)
  far_px <- rod(45, 45)
  # The same physical length occupies a different number of pixels, so no
  # single pixels-per-cm factor can be right for both.
  expect_gt(abs(near_px - far_px) / near_px, 0.1)

  # The homography recovers both correctly.
  check <- function(x0, y0) {
    a <- project_with(H, x0, y0)
    b <- project_with(H, x0 + 10, y0)
    ra <- rectify_coordinates(calibration, a$x, a$y)
    rb <- rectify_coordinates(calibration, b$x, b$y)
    sqrt((ra$x - rb$x)^2 + (ra$y - rb$y)^2)
  }
  expect_equal(check(5, 5), 10, tolerance = 1e-8)
  expect_equal(check(45, 45), 10, tolerance = 1e-8)
})

test_that("degenerate corners are rejected rather than silently fitted", {
  collinear <- data.frame(x = c(0, 1, 2, 3), y = c(0, 1, 2, 3))
  expect_error(arena_calibration(collinear, 10, 10))

  duplicated_corner <- data.frame(x = c(0, 0, 1, 1), y = c(0, 0, 1, 1))
  expect_error(arena_calibration(duplicated_corner, 10, 10))
})

test_that("calibration validates its arguments", {
  corners <- unit_square_corners()
  expect_error(arena_calibration(corners, -1, 10), "width")
  expect_error(arena_calibration(corners, 10, 0), "height")
  expect_error(arena_calibration(corners[1:3, ], 10, 10), "exactly four corners")
  expect_error(arena_calibration(corners, 10, 10, units = "furlong"), "Unsupported")
  expect_error(
    arena_calibration(data.frame(x = c(0, 1, 1, NA), y = c(1, 1, 0, 0)), 10, 10),
    "finite"
  )
})

test_that("rectify_tracking transforms landmarks and records provenance", {
  n <- 10
  corners <- data.frame(x = c(0, 500, 500, 0), y = c(500, 500, 0, 0))
  calibration <- arena_calibration(corners, 50, 50)

  points <- list(bodycentre = cbind(x = rep(250, n), y = rep(250, n)))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)

  rectified <- rectify_tracking(tracking, calibration)
  expect_equal(rectified$data$bodycentre$x, rep(25, n))
  expect_equal(rectified$data$bodycentre$y, rep(25, n))
  expect_equal(rectified$distance.units, "cm")
  expect_equal(rectified$arena_calibration$width, 50)
  expect_equal(rectified$median.data["bodycentre", "x"], 25)

  # Rectifying twice would treat centimetres as pixels.
  expect_error(rectify_tracking(rectified, calibration), "already calibrated")
})

test_that("arena violations are reported against the calibrated bounds", {
  n <- 6
  corners <- data.frame(x = c(0, 500, 500, 0), y = c(500, 500, 0, 0))
  calibration <- arena_calibration(corners, 50, 50)
  points <- list(
    bodycentre = cbind(
      x = c(250, 250, 250, 250, 5000, NA),
      y = c(250, 250, 250, 250, 250, NA)
    )
  )
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)
  rectified <- rectify_tracking(tracking, calibration)

  report <- arena_violation_report(rectified, "bodycentre", tolerance = 1)
  expect_equal(report$observed_frames, 5L)
  expect_equal(report$outside_arena_frames, 1L)
  expect_equal(report$outside_arena_fraction, 1 / 5)
})

test_that("points behind the camera plane become NA rather than reflections", {
  H <- matrix(c(1, 0, 0, 0, 1, 0, 0, 1, 1), nrow = 3, byrow = TRUE)
  mapped <- apply_homography(H, c(0, 0), c(0, -2))
  expect_equal(mapped$x[1], 0)
  expect_true(is.na(mapped$x[2]))
  expect_true(is.na(mapped$y[2]))
})
