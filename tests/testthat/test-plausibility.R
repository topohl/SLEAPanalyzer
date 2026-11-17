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
  expect_error(episodes_on_clip(transform(ep, end_frame = 118L), cut$clip_row, 30), "not aligned")
})

hand_events <- function() {
  data.frame(file = "synthetic", kind = c("start", "interior", "end"), first_frame = c(0, 40, 100),
             last_frame = c(2, 49, NA), evidence = "", reviewer = "", stringsAsFactors = FALSE)
}

test_that("hand events become leave-out windows by one rule", {
  windows <- hand_event_windows(hand_events(), fps = 30)
  # clip frame k is at k / fps: a start hand goes from frame 0, an end hand to the end of the recording
  expect_equal(windows$start_s, c(0, 40 / 30, 100 / 30))
  expect_equal(windows$end_s, c(2 / 30, 49 / 30, Inf))
  expect_equal(unique(windows$decision), "leave_out")
  expect_equal(unique(windows$target), "all")
  expect_match(windows$reason[3], "clip frames 100 to the end")
  # a start hand seen only from frame 1 still starts the phase once it has gone
  expect_equal(hand_event_windows(transform(hand_events()[1, ], first_frame = 1), 30)$start_s, 0)

  tracking <- climb_fixture()                            # 120 frames
  cut <- apply_leave_out(tracking, rbind(read_qc_decisions(NULL), windows), "synthetic")
  expect_equal(cut$clip_row, setdiff(1:120, c(1:3, 41:50, 101:120)))   # clip frame k is row k + 1
  expect_equal(cut$left_out_s, 33 / 30)
  expect_equal(cut$tracking$joins, 38L)                  # the interior hand leaves a join, the others none
  late <- hand_event_windows(transform(hand_events()[3, ], first_frame = 150), 30)
  expect_error(apply_leave_out(tracking, late, "synthetic"), "after the recording ends")
  # windows that leave no frame stop the run: an end hand at frame 0, a start hand past the end
  at_0 <- hand_event_windows(transform(hand_events()[3, ], first_frame = 0), 30)
  expect_error(apply_leave_out(tracking, at_0, "synthetic"), "leaves no frame")
  past <- hand_event_windows(transform(hand_events()[1, ], last_frame = 150), 30)
  expect_error(apply_leave_out(tracking, past, "synthetic"), "leaves no frame")

  path <- tempfile(fileext = ".csv")
  utils::write.csv(hand_events(), path, row.names = FALSE, na = "")
  expect_equal(hand_event_windows(read_hand_events(path), 30), windows)
  expect_equal(nrow(read_hand_events(file.path(tempdir(), "absent.csv"))), 0L)
  expect_equal(nrow(hand_event_windows(read_hand_events(NULL), 30)), 0L)
})

test_that("hand events are checked", {
  ev <- hand_events()[2, ]
  expect_error(hand_event_windows(transform(ev, kind = "middle"), 30), "Unknown hand event kind")
  expect_error(hand_event_windows(transform(ev, first_frame = 50), 30), "first_frame <= last_frame")
  expect_error(hand_event_windows(transform(ev, last_frame = NA), 30), "first_frame <= last_frame")
  expect_error(hand_event_windows(transform(ev, kind = "end"), 30), "no last_frame")
  expect_error(hand_event_windows(transform(ev, first_frame = 40.5), 30), "whole frame number")
  expect_error(hand_event_windows(transform(ev, first_frame = -1), 30), "whole frame number")
  expect_error(hand_event_windows(transform(ev, file = ""), 30), "a file in every row")
  expect_error(hand_event_windows(ev[, names(ev) != "reviewer"], 30), "lack column")
  expect_error(hand_event_windows(ev, 0), "fps")
  path <- tempfile(fileext = ".csv")
  utils::write.csv(transform(ev, kind = "end"), path, row.names = FALSE, na = "")
  expect_error(read_hand_events(path), "no last_frame")
})

test_that("a window inside a recording leaves a join that is no movement", {
  # 60 frames: still at 10 cm (1-25), walking 10 cm (26-35, left out), still at 20 cm (36-60)
  body <- rbind(xy(10, 12, 25), cbind(x = seq(11, 20, length.out = 10), y = rep(12, 10)), xy(20, 12, 25))
  nose <- cbind(x = body[, "x"] + 3, y = body[, "y"])
  tracking <- on_cup_tracking(body, nose, rep(1, 60))
  lo <- data.frame(file = "synthetic", target = "all", start_s = 25 / 30, end_s = 34 / 30,
                   decision = "leave_out", reason = "", reviewer = "", stringsAsFactors = FALSE)
  cut <- apply_leave_out(tracking, lo, "synthetic")
  expect_equal(cut$clip_row, c(1:25, 36:60))
  expect_equal(get_tracking_joins(cut$tracking), 26L)    # the first kept row after the window
  # the same rows as a recording of their own, where the 10 cm between rows 25 and 26 is a step
  own <- on_cup_tracking(body[cut$clip_row, ], nose[cut$clip_row, ], rep(1, 50))
  expect_equal(get_tracking_joins(own), integer())

  moved <- CalculateMovement(cut$tracking, movement_cutoff = 5, integration_period = 0)
  reference <- CalculateMovement(own, movement_cutoff = 5, integration_period = 0)
  expect_true(is.na(moved$data$bodycentre$speed[26]))
  after <- PointLocomotionReport(moved, "bodycentre"); before <- PointLocomotionReport(reference, "bodycentre")
  expect_equal(before$bodycentre.raw.distance, 10)
  expect_equal(after$bodycentre.raw.distance, 0)
  expect_equal(after$bodycentre.time.moving, 0)
  same <- c("bodycentre.total.time", "bodycentre.valid.time", "bodycentre.valid.fraction")
  expect_equal(after[same], before[same])
  columns <- function(t) unname(as.matrix(t$data$bodycentre[-26, c("x", "y", "speed")]))
  expect_equal(columns(moved), columns(reference))       # every other row and step as before
  # nor is it an implausible step
  expect_equal(displacement_qc(own, "bodycentre", max_speed = 100)$implausible_intervals, 1L)
  expect_equal(displacement_qc(cut$tracking, "bodycentre", max_speed = 100)$implausible_intervals, 0L)

  # no interpolation bridges the join: a frame missing on either side stays missing
  gap <- tracking
  gap$data$bodycentre[c(25, 36), c("x", "y")] <- NA
  bridged <- interpolate_tracking(apply_leave_out(gap, lo, "synthetic")$tracking, "bodycentre", max_gap_s = 0.2)
  expect_true(all(is.na(bridged$data$bodycentre$x[25:26])))
  own_gap <- own
  own_gap$data$bodycentre[25:26, c("x", "y")] <- NA
  expect_false(anyNA(interpolate_tracking(own_gap, "bodycentre", max_gap_s = 0.2)$data$bodycentre$x))
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

test_that("a reviewed not_hidden decision clears a hidden stretch and marks it", {
  cups <- c(left = "socl", right = "socr")
  tracking <- hidden_fixture()
  found <- detect_hidden(tracking)                       # clip frames 31-90 at the right cup, 1.00 to 3.00 s
  decision <- data.frame(file = "synthetic", target = "socr", start_s = 1, end_s = 3, decision = "not_hidden",
                         reason = "", reviewer = "", stringsAsFactors = FALSE)
  cleared <- apply_not_hidden(found, decision, "synthetic", cups, 30)
  expect_false(any(unlist(cleared$masks)))
  expect_equal(cleared$episodes$decision, "not_hidden")
  expect_identical(apply_not_hidden(found, decision, "other", cups, 30), found)
  # a stretch reaching well beyond the window is not the one reviewed
  expect_error(apply_not_hidden(found, transform(decision, end_s = 1.5), "synthetic", cups, 30), "review it again")
  # a window that meets no hidden stretch warns and clears nothing
  expect_warning(kept <- apply_not_hidden(found, transform(decision, target = "socl"), "synthetic", cups, 30),
                 "meets no hidden stretch")
  expect_equal(which(kept$masks$right), 31:90)
  expect_equal(kept$episodes$decision, "hidden")
  expect_error(apply_not_hidden(found, transform(decision, target = "objL"), "synthetic", cups, 30), "unknown target")

  # on the clip's clock after a leave-out at the start, longer than the 0.5 s an edge may move: counted on
  # the cut rows, the stretch would start 0.73 s before the window and the run would stop
  lo <- transform(decision, target = "all", start_s = 0, end_s = 0.7, decision = "leave_out")
  cut <- apply_leave_out(tracking, lo, "synthetic")
  expect_equal(cut$clip_row[1], 23L)
  on_cut <- apply_not_hidden(detect_hidden(cut$tracking), rbind(lo, decision), "synthetic", cups, 30, clip_row = cut$clip_row)
  expect_false(any(on_cut$masks$right))
  expect_equal(on_cut$episodes$decision, "not_hidden")
  expect_error(apply_not_hidden(detect_hidden(cut$tracking), rbind(lo, decision), "synthetic", cups, 30), "review it again")

  # not_hidden leaves the on-cup mask alone, and the on-cup decisions leave the hidden one alone
  expect_identical(apply_qc_decisions(found$masks$right, decision, "synthetic", "socr", 30), found$masks$right)
  not_on <- transform(decision, decision = "not_on_object")
  expect_identical(apply_not_hidden(found, not_on, "synthetic", cups, 30), found)

  path <- tempfile(fileext = ".csv")
  utils::write.csv(decision, path, row.names = FALSE)
  expect_equal(read_qc_decisions(path)$decision, "not_hidden")
  utils::write.csv(transform(decision, target = "all"), path, row.names = FALSE)
  expect_error(read_qc_decisions(path), "target 'all'")
})
