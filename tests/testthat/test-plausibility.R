# Frame-level plausibility checks: apparent size, on-object detection,
# reviewed decisions and the SocP contact mask. Synthetic trajectories in cm.

# A 44 x 24 cm arena with the right cup 4.5 cm from both walls of the
# bottom-right corner, as in the Exp9 social-preference arena.
on_cup_tracking <- function(body, nose, size, fps = 30) {
  n <- nrow(body)
  at <- function(x, y) cbind(x = rep(x, n), y = rep(y, n))
  # Every default size pair present (ear-ear, nose-ear, flank-flank) scales
  # with the apparent size.
  around <- function(origin, dx, dy) {
    cbind(x = origin[, "x"] + dx * size, y = origin[, "y"] + dy * size)
  }
  points <- list(
    tl = at(0, 24), tr = at(44, 24), br = at(44, 0), bl = at(0, 0),
    socl = at(4.5, 19.5), socr = at(39.5, 4.5),
    nose = nose, bodycentre = body,
    leftEar = around(nose, -0.5, -0.5), rightEar = around(nose, 0.5, -0.5),
    leftFlank = around(body, 0, -1), rightFlank = around(body, 0, 1)
  )
  make_tracking(frames = 0:(n - 1), fps = fps, points = points, distance_unit = "cm")
}

xy <- function(x, y, n) cbind(x = rep(x, n), y = rep(y, n))

# 120 frames: floor away (1-30), on top of the right cup and raised (31-60),
# still raised at the lid's front edge (61-75), floor contact (76-120).
climb_fixture <- function(raised = 1.2, on_frames = 30) {
  floor_away <- 30; edge <- 15; rest <- 120 - floor_away - on_frames - edge
  body <- rbind(xy(25, 12, floor_away), xy(38, 3, on_frames), xy(34, 6, edge), xy(34, 6, rest))
  nose <- rbind(xy(28, 12, floor_away), xy(37.5, 6.5, on_frames), xy(37, 6, edge), xy(37, 6, rest))
  size <- c(rep(1, floor_away), rep(raised, on_frames + edge), rep(1, rest))
  on_cup_tracking(body, nose, size)
}

detect_right <- function(tracking, ...) {
  detect_on_object(tracking, "socr", c("tl", "tr", "br", "bl"), footprint_radius = 4,
                   beyond = 2, max_distance = 10, min_apparent_size = 1.08,
                   merge_gap_s = 0.5, min_episode_s = 0.2, size_window_s = 0.5, ...)
}

test_that("apparent body size is relative to the recording median", {
  tracking <- climb_fixture()
  size <- apparent_body_size(tracking)
  expect_equal(size[1:30], rep(1, 30))
  expect_equal(size[31:75], rep(1.2, 45))
  expect_error(apparent_body_size(tracking, pairs = list(c("leftEar", "rightEar"))), "at least 2")
})

test_that("position relative to an object points toward the nearest corner", {
  tracking <- on_cup_tracking(rbind(xy(42.5, 1.5, 1), xy(34, 6, 1)), xy(30, 6, 2), c(1, 1))
  position <- object_relative_position(tracking, "socr", c("tl", "tr", "br", "bl"))
  expect_equal(position$distance, c(sqrt(18), sqrt(5.5^2 + 1.5^2)))
  expect_equal(position$beyond, c(sqrt(18), -(5.5 + 1.5) / sqrt(2)))
})

test_that("a raised animal over the cup is on the object, extended along the lid edge", {
  result <- detect_right(climb_fixture())
  expect_equal(which(result$on_object), 31:75)
  expect_equal(nrow(result$episodes), 1L)
  expect_equal(result$episodes$decision, "on_object")
  expect_equal(result$episodes$start_frame, 31L)
  expect_equal(result$episodes$end_frame, 75L)
  expect_equal(result$episodes$candidate_frames, 30L)
})

test_that("a normal-size animal in the same place is reported, not accepted", {
  result <- detect_right(climb_fixture(raised = 1))
  expect_false(any(result$on_object))
  expect_equal(result$episodes$decision, "normal_size")
})

test_that("an episode shorter than the minimum is reported, not accepted", {
  result <- detect_right(climb_fixture(on_frames = 3))
  expect_false(any(result$on_object))
  expect_equal(result$episodes$decision, "too_short")
})

test_that("reviewed decisions override the automatic mask", {
  tracking <- climb_fixture()
  mask <- detect_right(tracking)$on_object
  decisions <- data.frame(file = "synthetic", target = "socr", start_s = 1, end_s = 1.5,
                          decision = "not_on_object", reason = "", reviewer = "", stringsAsFactors = FALSE)
  forced <- apply_qc_decisions(mask, decisions, "synthetic", "socr", fps = 30)
  # Frame i lies at (i - 1) / 30 s: frames 31 to 46 span 1.0 to 1.5 s.
  expect_equal(which(forced), 47:75)
  expect_identical(apply_qc_decisions(mask, decisions, "other", "socr", 30), mask)

  path <- tempfile(fileext = ".csv")
  utils::write.csv(transform(decisions, decision = "maybe"), path, row.names = FALSE)
  expect_error(read_qc_decisions(path), "Unknown QC decision")
  expect_equal(nrow(read_qc_decisions(file.path(tempdir(), "absent.csv"))), 0L)
})

test_that("SocP does not count time on top of a cup as contact or proximity", {
  tracking <- climb_fixture()
  mask <- detect_right(tracking)$on_object
  before <- compute_socp_metrics(tracking, "L", fps = 30)
  after <- compute_socp_metrics(tracking, "L", fps = 30,
                                on_cup = list(left = rep(FALSE, 120), right = mask))
  # Contact with the right cup from frame 31 to 120; 31 to 75 are on top of it.
  expect_equal(before$summary$contactRight, 90 / 30)
  expect_equal(after$summary$contactRight, 45 / 30)
  expect_equal(after$summary$contactRemovedRight, 45 / 30)
  expect_equal(after$summary$onCupRight, 45 / 30)
  expect_equal(after$summary$onCupLeft, 0)
  # Metadata "L" maps the right cup to novel.
  expect_equal(after$summary$contactNovel, after$summary$contactRight)
  expect_equal(after$summary$onCupNovel, after$summary$onCupRight)
  # On-cup frames stay valid observation time.
  expect_equal(after$summary$validTime, before$summary$validTime)
  expect_equal(after$summary$latencyRight, 75 / 30)
  expect_error(compute_socp_metrics(tracking, "L", fps = 30, on_cup = list(right = mask)), "on_cup")
})

test_that("the report-only checks find jumps, distorted skeletons and points outside the arena", {
  n <- 60
  body <- xy(20, 12, n); body[31, ] <- c(20, 23); body[40:50, ] <- cbind(50, 12)
  tracking <- on_cup_tracking(body, cbind(x = body[, "x"] + 3, y = body[, "y"]), rep(1, n))
  jumps <- jump_episodes(tracking, "bodycentre", max_speed = 100)
  expect_true(all(c(31L, 32L, 40L) %in% unlist(Map(seq, jumps$start_frame, jumps$end_frame))))
  outside <- outside_arena_episodes(tracking, "bodycentre", c("tl", "tr", "br", "bl"))
  expect_equal(outside$start_frame, 40L)
  expect_equal(outside$end_frame, 50L)

  tracking$data$tailbase <- tracking$data$bodycentre
  tracking$data$tailbase$x <- tracking$data$tailbase$x - 3
  tracking$data$tailbase$x[10:20] <- tracking$data$tailbase$x[10:20] - 10
  skeleton <- skeleton_episodes(tracking, "nose", "tailbase", tolerance = 0.5)
  expect_equal(skeleton$start_frame, 10L)
  expect_equal(skeleton$end_frame, 20L)
})
