# Landmark repair: a misplaced landmark is found by the skeleton and the scored
# landmarks are rebuilt from rigid neighbours. A synthetic mouse in cm walks a
# circle with known body geometry.

synthetic_mouse <- function(n = 300, fps = 30, jitter = 0.03, seed = 1) {
  set.seed(seed)
  shape <- rbind(nose = c(4, 0), leftEar = c(2.5, 0.8), rightEar = c(2.5, -0.8), spine1 = c(1.5, 0),
                 leftArm = c(1.5, 1), rightArm = c(1.5, -1), bodycentre = c(0, 0), leftFlank = c(0, 1.4),
                 rightFlank = c(0, -1.4), spine2 = c(-1.5, 0), tailbase = c(-2.7, 0))
  t <- seq_len(n) / n * 2 * pi
  centre <- cbind(22 + 8 * cos(t), 12 + 8 * sin(t)); heading <- t + pi / 2
  points <- lapply(rownames(shape), function(p) {
    s <- shape[p, ]
    cbind(x = centre[, 1] + s[1] * cos(heading) - s[2] * sin(heading) + stats::rnorm(n, 0, jitter),
          y = centre[, 2] + s[1] * sin(heading) + s[2] * cos(heading) + stats::rnorm(n, 0, jitter))
  })
  names(points) <- rownames(shape)
  list(tracking = make_tracking(frames = 0:(n - 1), fps = fps, points = points, distance_unit = "cm"), truth = points)
}

test_that("a landmark far from the rest of the skeleton is flagged, and only that one", {
  m <- synthetic_mouse()
  tr <- m$tracking
  tr$data$nose$x[50:80] <- 2; tr$data$nose$y[50:80] <- 2          # the nose on something outside the animal
  flags <- skeleton_outliers(tr)
  expect_equal(which(flags[, "nose"]), 50:80)
  expect_false(any(flags[, setdiff(colnames(flags), "nose")]))
})

test_that("the misplaced nose is rebuilt from rigid neighbours, labelled reconstructed", {
  m <- synthetic_mouse()
  tr <- m$tracking
  tr$data$nose$x[50:80] <- 2; tr$data$nose$y[50:80] <- 2
  fixed <- repair_landmarks(tr)
  err <- sqrt((fixed$data$nose$x[50:80] - m$truth$nose[50:80, "x"])^2 + (fixed$data$nose$y[50:80] - m$truth$nose[50:80, "y"])^2)
  expect_lt(max(err), 0.25)
  expect_true(all(as.character(fixed$data$nose$status[50:80]) == "reconstructed"))
  expect_true(all(as.character(fixed$data$nose$status[-(50:80)]) == "observed"))
  expect_equal(fixed$data$leftEar$x, tr$data$leftEar$x)                 # nothing else touched
  s <- fixed$repair$summary
  expect_equal(s$flagged_s[s$landmark == "nose"], 31 / 30)
  expect_equal(s$reconstructed_s[s$landmark == "nose"], 31 / 30)
  expect_equal(s$still_missing_s[s$landmark == "nose"], 0)
  expect_lt(fixed$repair$models$nose$error_p95[1], 0.2)
  e <- fixed$repair$episodes
  expect_equal(e$start_frame[e$target == "nose"], 50L)
  expect_equal(e$end_frame[e$target == "nose"], 80L)
})

test_that("a collapsed skeleton is removed, not rebuilt", {
  m <- synthetic_mouse()
  tr <- m$tracking
  set.seed(9)
  for (p in names(tr$data)) { tr$data[[p]]$x[100:110] <- stats::runif(11, 0, 44); tr$data[[p]]$y[100:110] <- stats::runif(11, 0, 24) }
  fixed <- repair_landmarks(tr)
  expect_true(all(fixed$repair$collapse[100:110]))
  expect_false(any(fixed$repair$collapse[-(100:110)]))
  expect_true(all(is.na(fixed$data$nose$x[100:110])))
  expect_true(all(as.character(fixed$data$bodycentre$status[100:110]) == "invalid"))
  expect_equal(fixed$repair$summary$collapse_s[1], 11 / 30)
  expect_true("skeleton_collapse" %in% fixed$repair$episodes$check)
})

test_that("a split skeleton keeps its consistent majority; without one it collapses", {
  m <- synthetic_mouse()
  tr <- m$tracking
  head <- c("nose", "leftEar", "rightEar", "spine1", "rightArm")       # 5 of 11 moved together, off the animal
  for (p in head) { tr$data[[p]]$x[100:130] <- tr$data[[p]]$x[100:130] - 15; tr$data[[p]]$y[100:130] <- tr$data[[p]]$y[100:130] + 40 }
  expect_true(all(skeleton_outliers(tr)[100:130, ]))                   # the vote alone flags every landmark
  fixed <- repair_landmarks(tr)
  expect_false(any(fixed$repair$collapse))
  expect_true(all(fixed$repair$split[100:130]))
  expect_false(any(fixed$repair$split[-(100:130)]))
  expect_true(all(is.na(fixed$data$leftEar$x[100:130])))
  expect_equal(fixed$data$tailbase$x, tr$data$tailbase$x)
  err <- sqrt((fixed$data$nose$x[100:130] - m$truth$nose[100:130, "x"])^2 + (fixed$data$nose$y[100:130] - m$truth$nose[100:130, "y"])^2)
  expect_lt(max(err), 0.25)
  expect_true(all(as.character(fixed$data$nose$status[100:130]) == "reconstructed"))
  expect_true("skeleton_split" %in% fixed$repair$episodes$check)
  tr$data$tailbase$x[100:130] <- 40; tr$data$tailbase$y[100:130] <- 2  # now 5 against 5 and one: no majority
  even <- repair_landmarks(tr)
  expect_true(all(even$repair$collapse[100:130]))
  expect_false(any(even$repair$split))
})

test_that("landmarks far outside the arena are removed before the skeleton vote", {
  m <- synthetic_mouse()
  n <- nrow(m$truth$nose)
  corner <- function(x, y) cbind(x = rep(x, n), y = rep(y, n))
  points <- c(m$truth, list(tl = corner(0, 24), tr = corner(44, 24), br = corner(44, 0), bl = corner(0, 0)))
  out <- c("nose", "leftEar", "rightEar", "spine1", "rightArm", "leftArm")   # 6 of 11 stuck on one spot 30 cm outside
  for (p in out) points[[p]][100:130, ] <- cbind(-30 + 0.3 * match(p, out), 10)
  tr <- make_tracking(frames = 0:(n - 1), fps = 30, points = points, distance_unit = "cm")
  expect_equal(distance_outside_polygon(c(22, -3, 50, NA), c(12, 12, 30, 1), rbind(c(0, 24), c(44, 24), c(44, 0), c(0, 0))), c(0, 3, sqrt(72), NA))
  fixed <- repair_landmarks(tr, corners = c("tl", "tr", "br", "bl"), max_outside = 15)
  expect_false(any(fixed$repair$collapse))
  err <- sqrt((fixed$data$nose$x[100:130] - m$truth$nose[100:130, "x"])^2 + (fixed$data$nose$y[100:130] - m$truth$nose[100:130, "y"])^2)
  expect_lt(max(err), 0.25)
  expect_equal(fixed$data$tailbase$x, tr$data$tailbase$x)
  s <- fixed$repair$summary
  expect_equal(s$outside_s[s$landmark == "leftEar"], 31 / 30)
  expect_equal(s$outside_s[s$landmark == "tailbase"], 0)
  expect_true(all(fixed$repair$frames_flagged[100:130]))
  expect_true("outside_arena" %in% fixed$repair$episodes$check)
  without <- repair_landmarks(tr)                                     # the vote alone keeps the larger, wrong group
  expect_true(all(is.na(without$data$tailbase$x[100:130])))
  for (p in names(m$truth)) points[[p]][200:205, ] <- cbind(60 + 0.3 * match(p, names(m$truth)), 50)
  gone <- repair_landmarks(make_tracking(frames = 0:(n - 1), fps = 30, points = points, distance_unit = "cm"),
                           corners = c("tl", "tr", "br", "bl"), max_outside = 15)
  expect_true(all(gone$repair$collapse[200:205]))                     # every landmark outside: a collapsed skeleton
  expect_true(all(is.na(gone$data$nose$x[200:205])))
})

test_that("rebuilt frames survive interpolation and count as valid but not as observed", {
  m <- synthetic_mouse()
  tr <- m$tracking
  tr$data$nose$x[50:80] <- 2; tr$data$nose$y[50:80] <- 2
  tr$data$nose$x[200:203] <- NA; tr$data$nose$y[200:203] <- NA          # untracked, but anchors present: rebuilt too
  tr$data$bodycentre$x[150:152] <- NA; tr$data$bodycentre$y[150:152] <- NA
  fixed <- interpolate_tracking(repair_landmarks(tr), landmarks = c("nose", "bodycentre"), max_gap_s = 0.2)
  expect_true(all(as.character(fixed$data$nose$status[c(50:80, 200:203)]) == "reconstructed"))
  expect_true(all(as.character(fixed$data$bodycentre$status[150:152]) == "reconstructed"))
  expect_true(all(landmark_validity(fixed, c("nose", "bodycentre"))))
  expect_equal(sum(!landmark_validity(fixed, c("nose", "bodycentre"), reconstructed_is_valid = FALSE)), 31 + 4 + 3)
  qc <- tracking_qc_report(fixed, landmarks = c("nose", "bodycentre"))
  expect_equal(qc$reconstructed_time_s, 38 / 30)
  expect_equal(qc$interpolated_time_s, 0)
  expect_equal(qc$observed_time_s, (300 - 38) / 30)
  expect_equal(qc$valid_time_s, 10)
})

test_that("anchor pairs with too large an error are not used", {
  m <- synthetic_mouse(jitter = 0.6)
  tr <- m$tracking
  tr$data$nose$x[50:80] <- 2; tr$data$nose$y[50:80] <- 2
  fixed <- repair_landmarks(tr, max_anchor_error = 0.1)
  expect_true(all(is.na(fixed$data$nose$x[50:80])))
  expect_equal(fixed$repair$summary$reconstructed_s[1], 0)
})
