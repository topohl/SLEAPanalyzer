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
  expect_equal(which(mask), 31:75)                       # frame i lies at (i - 1) / 30 s: 1.00 to 2.47 s
  decisions <- data.frame(file = "synthetic", target = "socr", start_s = 1.1, end_s = 2.3,
                          decision = "not_on_object", reason = "", reviewer = "", stringsAsFactors = FALSE)
  # the reviewed episode, its edges moved by up to 0.5 s: cleared as a whole
  expect_false(any(apply_qc_decisions(mask, decisions, "synthetic", "socr", fps = 30)))
  expect_identical(apply_qc_decisions(mask, decisions, "other", "socr", 30), mask)
  # an episode reaching well beyond the reviewed window is not the one reviewed
  expect_error(apply_qc_decisions(mask, transform(decisions, end_s = 1.5), "synthetic", "socr", 30), "review it again")
  forced <- apply_qc_decisions(mask, transform(decisions, start_s = 3, end_s = 3.2, decision = "on_object"), "synthetic", "socr", 30)
  expect_equal(which(forced), c(31:75, 91:97))           # 3.0 to 3.2 s forced on, the detected episode kept

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

test_that("leave-out windows are read, need target 'all' and remove their frames", {
  decisions <- data.frame(file = "synthetic", target = "all", start_s = 1, end_s = 2,
                          decision = "leave_out", reason = "", reviewer = "", stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".csv")
  utils::write.csv(decisions, path, row.names = FALSE)
  expect_equal(read_qc_decisions(path)$decision, "leave_out")
  utils::write.csv(transform(decisions, target = "socl"), path, row.names = FALSE)
  expect_error(read_qc_decisions(path), "target 'all'")
  utils::write.csv(transform(decisions, decision = "on_object"), path, row.names = FALSE)
  expect_error(read_qc_decisions(path), "target 'all'")

  tracking <- climb_fixture()                            # 120 frames: frame i at (i - 1) / 30 s
  cut <- apply_leave_out(tracking, decisions, "synthetic")
  expect_equal(cut$left_out_s, 31 / 30)                  # frames 31 to 61, 1.00 to 2.00 s
  expect_equal(length(cut$tracking$frames), 89L)
  expect_equal(cut$tracking$frames, tracking$frames[-(31:61)])
  expect_equal(nrow(cut$tracking$data$nose), 89L)
  expect_equal(cut$clip_row, setdiff(1:120, 31:61))       # each kept row's row in the clip
  expect_identical(apply_leave_out(tracking, decisions, "other")$tracking, tracking)
  expect_equal(apply_leave_out(tracking, decisions, "other")$clip_row, 1:120)
  # on-object decisions are not leave-outs
  expect_identical(apply_leave_out(tracking, transform(decisions, decision = "on_object"), "synthetic")$tracking, tracking)
  expect_error(apply_leave_out(tracking, transform(decisions, start_s = 5, end_s = 6), "synthetic"), "after the recording ends")
})

test_that("after a leave-out, decisions and episode tables stay on the clip's clock", {
  tracking <- climb_fixture()
  mask <- detect_right(tracking)$on_object               # clip frames 31-75
  lo <- data.frame(file = "synthetic", target = "all", start_s = 0, end_s = 0.07,
                   decision = "leave_out", reason = "", reviewer = "", stringsAsFactors = FALSE)
  cut <- apply_leave_out(tracking, lo, "synthetic")      # the first three frames go
  expect_equal(cut$clip_row, 4:120)
  kept <- mask[cut$clip_row]
  on <- transform(lo, target = "socr", start_s = 3, end_s = 3.2, decision = "on_object")
  forced <- apply_qc_decisions(kept, on, "synthetic", "socr", 30, clip_row = cut$clip_row)
  expect_equal(cut$clip_row[forced], c(31:75, 91:97))   # the same clip frames as without the leave-out
  # counted on the cut rows, the window would land three frames late
  expect_equal(cut$clip_row[apply_qc_decisions(kept, on, "synthetic", "socr", 30)], c(31:75, 94:100))
  expect_error(apply_qc_decisions(kept, on, "synthetic", "socr", 30, clip_row = 1:120), "not aligned")

  ep <- data.frame(start_frame = 28L, end_frame = 72L, start_s = 27 / 30, end_s = 72 / 30, duration_s = 45 / 30)
  moved <- episodes_on_clip(ep, cut$clip_row, 30)
  expect_equal(c(moved$start_frame, moved$end_frame), c(31L, 75L))
  expect_equal(c(moved$start_s, moved$end_s, moved$duration_s), c(1, 2.5, 1.5))
  expect_identical(episodes_on_clip(ep, 1:120, 30), ep)  # no leave-out: unchanged
  expect_identical(episodes_on_clip(ep[0, ], cut$clip_row, 30), ep[0, ])
})

# 120 frames near the right cup; frames 31-90 are a gap fill of the body (a
# straight line at constant speed between the observed frames 30 and 91).
hidden_fixture <- function(nose_in_fill = c("filled", "seen"), body_at = c(35, 6)) {
  n <- 120
  wobble <- 0.3 * sin(seq_len(n) * 1.3)
  body <- cbind(x = body_at[1] + wobble, y = body_at[2] + wobble / 2)
  nose <- cbind(x = body[, "x"] - 3, y = body[, "y"] + 1)
  fill <- 31:90
  line <- function(m) {
    for (j in 1:2) m[fill, j] <- m[30, j] + (m[91, j] - m[30, j]) * (fill - 30) / 61
    m
  }
  body <- line(body)
  if (match.arg(nose_in_fill) == "filled") nose <- line(nose)
  tracking <- on_cup_tracking(body, nose, rep(1, n))
  tracking$px.to.cm <- 0.05
  tracking
}

detect_hidden <- function(tracking, ...) {
  detect_hidden_at_objects(tracking, c(left = "socl", right = "socr"), min_s = 0.2, max_distance = 10,
                           nose_off_body = 10, ...)
}

test_that("a mouse lost beside a cup is hidden there; a seen nose or a far body is not", {
  found <- detect_hidden(hidden_fixture())
  expect_equal(which(found$masks$right), 31:90)
  expect_false(any(found$masks$left))
  expect_equal(found$episodes$decision, "hidden")
  expect_equal(found$episodes$target, "socr")
  expect_equal(found$episodes$start_frame, 31L)
  expect_equal(found$episodes$duration_s, 2)

  seen <- detect_hidden(hidden_fixture("seen"))
  expect_false(any(seen$masks$right))
  expect_equal(seen$episodes$decision, "nose_seen")

  away <- detect_hidden(hidden_fixture(body_at = c(22, 12)))
  expect_false(any(unlist(away$masks)))
  expect_equal(away$episodes$decision, "away")

  # a nose thrown far from the body is lost too
  thrown <- hidden_fixture("seen")
  thrown$data$nose$x[31:90] <- 1; thrown$data$nose$y[31:90] <- 23
  expect_equal(which(detect_hidden(thrown)$masks$right), 31:90)

  # a fill shorter than min_s is not reported
  expect_equal(nrow(detect_hidden_at_objects(hidden_fixture(), c(left = "socl", right = "socr"), min_s = 2.5,
                                             max_distance = 10, nose_off_body = 10)$episodes), 0L)
})

test_that("SocP counts hidden time as contact with that cup only, and as valid time", {
  tracking <- hidden_fixture()
  hidden <- detect_hidden(tracking)$masks
  # the nose is missing while hidden: without the rule those frames are invalid
  tracking$data$nose$x[31:90] <- NA; tracking$data$nose$y[31:90] <- NA
  before <- compute_socp_metrics(tracking, "L", fps = 30)
  after <- compute_socp_metrics(tracking, "L", fps = 30, hidden = hidden,
                                on_cup = list(left = rep(FALSE, 120), right = seq_len(120) %in% 41:50))
  expect_equal(after$summary$contactRight, before$summary$contactRight + 60 / 30)
  expect_equal(after$summary$hiddenRight, 2)
  expect_equal(after$summary$hiddenNovel, 2)             # "L": the right cup is novel
  expect_equal(after$summary$hiddenLeft, 0)
  expect_equal(after$summary$onCupRight, 0)              # hidden is never on the cup
  expect_equal(after$summary$validTime, 4)
  expect_equal(before$summary$validTime, 2)
  expect_false(any(after$events$proximity_right[31:90]))
  expect_false(any(after$events$contact_left[31:90]))
  expect_error(compute_socp_metrics(tracking, "L", fps = 30, hidden = list(right = hidden$right)), "hidden")
})
