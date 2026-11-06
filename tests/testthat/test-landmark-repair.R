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
  # the head went with the removed group: without ears the nose cannot be checked against the head, so it stays missing
  expect_true(all(is.na(fixed$data$nose$x[100:130])))
  expect_equal(fixed$repair$summary$refused_s[fixed$repair$summary$landmark == "nose"], 31 / 30)
  expect_true(all(is.finite(fixed$data$bodycentre$x[100:130])))
  without <- repair_landmarks(tr, head_check = FALSE)                  # without the check it is rebuilt from the body
  err <- sqrt((without$data$nose$x[100:130] - m$truth$nose[100:130, "x"])^2 + (without$data$nose$y[100:130] - m$truth$nose[100:130, "y"])^2)
  expect_lt(max(err), 0.25)
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
  fixed <- repair_landmarks(tr, corners = c("tl", "tr", "br", "bl"), max_outside = 15, head_check = FALSE)
  expect_false(any(fixed$repair$collapse))
  err <- sqrt((fixed$data$nose$x[100:130] - m$truth$nose[100:130, "x"])^2 + (fixed$data$nose$y[100:130] - m$truth$nose[100:130, "y"])^2)
  expect_lt(max(err), 0.25)
  expect_true(all(is.na(repair_landmarks(tr, corners = c("tl", "tr", "br", "bl"), max_outside = 15)$data$nose$x[100:130])))   # ears gone: not rebuilt
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

test_that("a nose short of a hidden head is put back at the median head length, along the head", {
  m <- synthetic_mouse()
  tr <- m$tracking
  tr$data$nose$status <- factor(rep("observed", 300), levels = FRAME_STATUS_LEVELS)
  ex <- (tr$data$leftEar$x + tr$data$rightEar$x) / 2; ey <- (tr$data$leftEar$y + tr$data$rightEar$y) / 2
  k <- 50:80                                                             # head under a lid: the nose seen at 40% of the head
  tr$data$nose$x[k] <- ex[k] + 0.4 * (m$truth$nose[k, "x"] - ex[k]); tr$data$nose$y[k] <- ey[k] + 0.4 * (m$truth$nose[k, "y"] - ey[k])
  out <- restore_head_length(tr)
  expect_equal(out$head_length$median, 1.5, tolerance = 0.02)
  err <- sqrt((out$data$nose$x[k] - m$truth$nose[k, "x"])^2 + (out$data$nose$y[k] - m$truth$nose[k, "y"])^2)
  expect_lt(max(err), 0.1)
  expect_equal(which(out$head_length$moved), k)                         # only the clearly short heads
  expect_equal(out$data$nose$x[-k], tr$data$nose$x[-k])
  expect_true(all(as.character(out$data$nose$status[k]) == "reconstructed"))
  expect_true(all(as.character(out$data$nose$status[-k]) == "observed"))
  expect_equal(out$data$leftEar$x, tr$data$leftEar$x)
})

test_that("a head near its median length is left as it is; a short head turned sideways is restored along its own direction", {
  m <- synthetic_mouse()
  tr <- m$tracking
  ex <- (tr$data$leftEar$x + tr$data$rightEar$x) / 2; ey <- (tr$data$leftEar$y + tr$data$rightEar$y) / 2
  bx <- ex - tr$data$bodycentre$x; by <- ey - tr$data$bodycentre$y; bl <- sqrt(bx^2 + by^2)
  side <- 60:70                                                          # short and at right angles to the body axis
  tr$data$nose$x[side] <- ex[side] - 0.6 * by[side] / bl[side]; tr$data$nose$y[side] <- ey[side] + 0.6 * bx[side] / bl[side]
  near <- 100:110                                                        # 0.9 of the median: not clearly short
  tr$data$nose$x[near] <- ex[near] + 0.9 * (m$truth$nose[near, "x"] - ex[near]); tr$data$nose$y[near] <- ey[near] + 0.9 * (m$truth$nose[near, "y"] - ey[near])
  out <- restore_head_length(tr)
  expect_false(any(out$head_length$moved[near]))
  expect_equal(out$data$nose$x[near], tr$data$nose$x[near])
  expect_true(all(out$head_length$moved[side]))
  turned <- sqrt((out$data$nose$x[side] - ex[side])^2 + (out$data$nose$y[side] - ey[side])^2)
  expect_equal(turned, rep(out$head_length$median, length(side)), tolerance = 1e-9)   # full length, still at right angles
  expect_lte(max(out$head_length$shift), out$head_length$median)
})

test_that("head length restoration leaves frames without both ears, and takes its median from observed noses only", {
  m <- synthetic_mouse()
  tr <- m$tracking
  tr$data$leftEar$x[10:20] <- NA; tr$data$leftEar$y[10:20] <- NA
  tr$data$nose$x[10:20] <- tr$data$rightEar$x[10:20]; tr$data$nose$y[10:20] <- tr$data$rightEar$y[10:20]
  tr$data$nose$status <- factor(rep("observed", 300), levels = FRAME_STATUS_LEVELS)
  ex <- (tr$data$leftEar$x + tr$data$rightEar$x) / 2; ey <- (tr$data$leftEar$y + tr$data$rightEar$y) / 2
  k <- 100:250                                                           # long stretch of short heads, not observed: no effect on the median
  tr$data$nose$x[k] <- ex[k] + 0.5 * (tr$data$nose$x[k] - ex[k]); tr$data$nose$y[k] <- ey[k] + 0.5 * (tr$data$nose$y[k] - ey[k])
  tr$data$nose$status[k] <- "reconstructed"
  out <- restore_head_length(tr)
  expect_equal(out$data$nose$x[10:20], tr$data$nose$x[10:20])
  expect_equal(out$head_length$median, 1.5, tolerance = 0.02)
  expect_true(all(out$head_length$moved[k]))
})

test_that("without an observed nose and both ears, head length restoration warns and changes nothing", {
  m <- synthetic_mouse()
  tr <- m$tracking
  tr$data$nose$status <- factor(rep("reconstructed", 300), levels = FRAME_STATUS_LEVELS)
  expect_warning(out <- restore_head_length(tr), "no frame with an observed nose")
  expect_equal(out$data$nose$x, tr$data$nose$x)
  expect_true(is.na(out$head_length$median))
  expect_false(any(out$head_length$moved))
})

test_that("a rebuilt nose that would contradict the head is refused, not placed", {
  m <- synthetic_mouse()
  tr <- m$tracking
  tr$data$nose$x[50:80] <- NA; tr$data$nose$y[50:80] <- NA
  ok <- repair_landmarks(tr)                                            # ears present and right: rebuilt in front of the head
  err <- sqrt((ok$data$nose$x[50:80] - m$truth$nose[50:80, "x"])^2 + (ok$data$nose$y[50:80] - m$truth$nose[50:80, "y"])^2)
  expect_lt(max(err), 0.25)
  sw <- tr; sw$data$leftEar <- tr$data$rightEar; sw$data$rightEar <- tr$data$leftEar   # ears swapped: the ear pair points backwards
  sw$data$leftEar$x[-(50:80)] <- tr$data$leftEar$x[-(50:80)]; sw$data$leftEar$y[-(50:80)] <- tr$data$leftEar$y[-(50:80)]
  sw$data$rightEar$x[-(50:80)] <- tr$data$rightEar$x[-(50:80)]; sw$data$rightEar$y[-(50:80)] <- tr$data$rightEar$y[-(50:80)]
  bad <- repair_landmarks(sw)
  ex <- (sw$data$leftEar$x + sw$data$rightEar$x) / 2; ey <- (sw$data$leftEar$y + sw$data$rightEar$y) / 2
  ax <- ex - sw$data$bodycentre$x; ay <- ey - sw$data$bodycentre$y
  hx <- bad$data$nose$x - ex; hy <- bad$data$nose$y - ey
  front <- (hx * ax + hy * ay)[50:80]
  expect_true(all(is.na(front) | front > 0))                           # whatever is placed lies in front of the head
  expect_gt(sum(is.finite(bad$data$nose$x[50:80])), 0)                 # the next anchor pair places it
})

test_that("a rebuilt nose keeps the head pose of the moment, fading to the typical pose deep in a gap", {
  m <- synthetic_mouse(n = 600, jitter = 0.01)
  tr <- m$tracking
  ex <- (tr$data$leftEar$x + tr$data$rightEar$x) / 2; ey <- (tr$data$leftEar$y + tr$data$rightEar$y) / 2
  th <- 0.6 * sin(2 * pi * seq_len(600) / 300)                         # the head turns slowly to and fro
  hx <- tr$data$nose$x - ex; hy <- tr$data$nose$y - ey
  truth <- cbind(ex + cos(th) * hx - sin(th) * hy, ey + sin(th) * hx + cos(th) * hy)
  tr$data$nose$x <- truth[, 1]; tr$data$nose$y <- truth[, 2]
  k <- 70:79                                                            # nose hidden for a third of a second, head turned
  tr$data$nose$x[k] <- NA; tr$data$nose$y[k] <- NA
  err <- function(x) sqrt((x$data$nose$x[k] - truth[k, 1])^2 + (x$data$nose$y[k] - truth[k, 2])^2)
  pose <- repair_landmarks(tr, pose_decay_s = 1); typical <- repair_landmarks(tr, pose_decay_s = 0)
  expect_true(all(as.character(pose$data$nose$status[k]) == "reconstructed"))
  expect_lt(max(err(pose)), 0.5 * stats::median(err(typical)))         # the pose of the moment, far closer than the typical pose
  expect_gt(stats::median(err(typical)), 0.4)                           # the typical pose misses the turn
  long <- 200:499                                                       # ten seconds: the middle falls back to the typical pose
  tr2 <- tr; tr2$data$nose$x[long] <- NA; tr2$data$nose$y[long] <- NA
  a <- repair_landmarks(tr2, pose_decay_s = 1); b <- repair_landmarks(tr2, pose_decay_s = 0)
  mid <- 330:370
  expect_lt(max(abs(a$data$nose$x[mid] - b$data$nose$x[mid])), 0.05)
})

test_that("the pose of the moment works at the start and end of a recording and for the body centre", {
  m <- synthetic_mouse()
  tr <- m$tracking
  ends <- c(1:10, 291:300)
  tr$data$nose$x[ends] <- NA; tr$data$nose$y[ends] <- NA
  tr$data$bodycentre$x[ends] <- NA; tr$data$bodycentre$y[ends] <- NA
  out <- repair_landmarks(tr, pose_decay_s = 1)
  expect_true(all(is.finite(out$data$nose$x[ends])))
  expect_true(all(as.character(out$data$nose$status[ends]) == "reconstructed"))
  expect_true(all(is.finite(out$data$bodycentre$x[ends])))
  s <- out$repair$summary
  expect_equal(s$refused_s[s$landmark == "nose"], 0)
  expect_equal(s$still_missing_s[s$landmark == "bodycentre"], 0)
  err <- sqrt((out$data$nose$x[ends] - m$truth$nose[ends, "x"])^2 + (out$data$nose$y[ends] - m$truth$nose[ends, "y"])^2)
  expect_lt(max(err), 0.25)
})

test_that("an anchor pair squeezed together at the edge of a gap does not carry its pose in", {
  m <- synthetic_mouse()
  tr <- m$tracking
  tr$data$rightEar$x[100] <- tr$data$leftEar$x[100] + 0.01; tr$data$rightEar$y[100] <- tr$data$leftEar$y[100]
  tr$data$nose$x[101:110] <- NA; tr$data$nose$y[101:110] <- NA
  out <- repair_landmarks(tr, pose_decay_s = 1)
  err <- sqrt((out$data$nose$x[101:110] - m$truth$nose[101:110, "x"])^2 + (out$data$nose$y[101:110] - m$truth$nose[101:110, "y"])^2)
  expect_true(all(is.finite(err)))
  expect_lt(max(err), 0.25)
})

test_that("without the pose of the moment the nose sits at the anchor pair's median position", {
  m <- synthetic_mouse()
  tr <- m$tracking
  tr$data$nose$x[50:80] <- NA; tr$data$nose$y[50:80] <- NA
  out <- repair_landmarks(tr, pose_decay_s = 0)
  best <- out$repair$models$nose[1, ]
  L <- list(leftEar = cbind(tr$data$leftEar$x, tr$data$leftEar$y), rightEar = cbind(tr$data$rightEar$x, tr$data$rightEar$y))
  L[[best$anchor_a]] <- cbind(tr$data[[best$anchor_a]]$x, tr$data[[best$anchor_a]]$y); L[[best$anchor_b]] <- cbind(tr$data[[best$anchor_b]]$x, tr$data[[best$anchor_b]]$y)
  expect_equal(unname(cbind(out$data$nose$x[50:80], out$data$nose$y[50:80])), unname(predict_from_anchor(L, best, 50:80)))
  expect_error(repair_landmarks(tr, pose_decay_s = -1))
  expect_error(repair_landmarks(tr, pose_decay_s = NA))
})
