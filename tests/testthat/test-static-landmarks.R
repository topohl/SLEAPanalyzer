test_that("static landmarks are held at their median over the frames in which they are apart", {
  n <- 10
  objL <- cbind(x = c(0, 0.2, -0.1, 0, 20, 20, 0.1, 0, NA, 0), y = rep(0, n))   # frames 5-6: the keypoint jumped onto objR
  objR <- cbind(x = rep(20, n), y = c(0, 0, 0.3, 0, 0, 0, 0, -0.2, 0, 0))
  tr <- make_tracking(frames = 0:(n - 1), points = list(objL = objL, objR = objR))
  out <- hold_static_landmarks(tr, c("objL", "objR"), min_separation = 10)
  apart <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE)
  expect_equal(out$data$objL$x, rep(median(objL[apart, "x"]), n))
  expect_equal(out$data$objR$y, rep(median(objR[apart, "y"]), n))
  expect_true(all(is.finite(out$data$objL$x)))                                  # missing frames get the fixed position too
  expect_equal(attr(out, "static_landmarks")$share_of_frames_used, mean(apart))
  expect_equal(out$median.data["objL", "x"], out$data$objL$x[1])
})

test_that("static landmarks refuse landmarks that are never apart, and a wrong number of landmarks", {
  tr <- make_tracking(frames = 0:4, points = list(objL = cbind(x = rep(1, 5), y = 0), objR = cbind(x = rep(2, 5), y = 0)))
  expect_error(hold_static_landmarks(tr, c("objL", "objR"), min_separation = 10), "never")
  expect_error(hold_static_landmarks(tr, "objL", min_separation = 10), "exactly two")
})
