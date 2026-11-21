# End-to-end run of the SocP batch script against synthetic data.

make_socp_fixture <- function(root, n = 300, fps = 30, contact_frames = 90,
                              climb_frames = 0, hidden_frames = 0, extra_lines = character()) {
  # 440 x 240 px arena for a 44 x 24 cm apparatus, so 10 px per cm.
  points <- list(
    tl = cbind(x = rep(0, n), y = rep(240, n)),
    tr = cbind(x = rep(440, n), y = rep(240, n)),
    br = cbind(x = rep(440, n), y = rep(0, n)),
    bl = cbind(x = rep(0, n), y = rep(0, n)),
    socl = cbind(x = rep(50, n), y = rep(120, n)),
    socr = cbind(x = rep(390, n), y = rep(120, n))
  )
  # The nose sits on the left stimulus for contact_frames, then in the middle.
  nose_x <- c(rep(50, contact_frames), rep(220, n - contact_frames))
  nose_y <- rep(120, n)
  points$nose <- cbind(x = nose_x, y = nose_y)
  points$bodycentre <- cbind(x = nose_x + 30, y = nose_y)

  on_cup_lines <- character()
  if (climb_frames > 0) {
    # Floor contact with the body 5 cm from the cup point, then the animal sits
    # on top of the left cup, its body 2 cm from the cup point and imaged 1.25
    # times larger, nose still within the contact distance.
    points$bodycentre[seq_len(contact_frames), "x"] <- 100
    climb <- contact_frames + seq_len(climb_frames)
    points$nose[climb, "x"] <- 50
    points$bodycentre[climb, "x"] <- 70
    size <- rep(1, n); size[climb] <- 1.25
    ears <- function(dx) cbind(x = points$nose[, "x"] + dx * 5 * size, y = points$nose[, "y"] - 5 * size)
    flank <- function(dy) cbind(x = points$bodycentre[, "x"], y = points$bodycentre[, "y"] + dy * 10 * size)
    points$leftEar <- ears(-1); points$rightEar <- ears(1)
    points$leftFlank <- flank(-1); points$rightFlank <- flank(1)
    on_cup_lines <- c(
      "on_cup_exclusion: true", "on_cup_footprint_radius_cm: 4", "on_cup_beyond_cm: 2",
      "on_cup_max_distance_cm: 10", "on_cup_min_apparent_size: 1.08",
      "on_cup_merge_gap_s: 0.5", "on_cup_min_episode_s: 0.2", "on_cup_size_window_s: 0.5"
    )
  }

  if (hidden_frames > 0) {
    # Tracked frames move a little, as real tracking does. The mouse comes to
    # the right cup (body 6 cm from the cup point, nose in contact) for 10
    # frames, is lost for hidden_frames (the export fills the body by a straight
    # line and the pose model puts the nose at the image corner), and is seen
    # again for 10 frames.
    wobble <- 0.5 * sin(seq_len(n) * 1.3)
    points$nose[, "x"] <- points$nose[, "x"] + wobble
    points$bodycentre[, "x"] <- points$bodycentre[, "x"] + wobble
    a <- 141; b <- a + 10 + hidden_frames + 9
    visit <- a:b
    points$bodycentre[visit, "x"] <- 330 + wobble[visit]
    points$nose[visit, "x"] <- 360 + wobble[visit]
    fill <- (a + 10):(a + 9 + hidden_frames)
    first <- min(fill) - 1; last <- max(fill) + 1
    for (j in c("x", "y")) {
      points$bodycentre[fill, j] <- points$bodycentre[first, j] +
        (points$bodycentre[last, j] - points$bodycentre[first, j]) * (fill - first) / (last - first)
    }
    points$nose[fill, "x"] <- 10 + wobble[fill]
    points$nose[fill, "y"] <- 10
  }

  write_tracking_fixture(
    file.path(root, "formatted", "HAB", "B001_socp.csv"), points, n
  )
  writeLines(c("Code\tNovelLoc", "B001\tR"), file.path(root, "novelLocHAB.txt"))
  writeLines(c("Code ID", "B001 mouse-2"), file.path(root, "animalIDCode.txt"))

  config_path <- file.path(root, "socp.yaml")
  writeLines(c(
    paste0("input_dir: ", file.path(root, "formatted")),
    paste0("output_dir: ", file.path(root, "output")),
    paste0("metadata_dir: ", root),
    paste0("animal_id_code_file: ", file.path(root, "animalIDCode.txt")),
    "novel_location_file_prefix: novelLoc",
    "phases: [HAB]",
    paste0("fps: ", fps),
    "arena_width_cm: 44",
    "arena_height_cm: 24",
    "arena_corner_names: [tl, tr, br, bl]",
    "max_interpolation_gap_s: 0.2",
    "likelihood_cutoff: ~",
    "max_plausible_speed_cm_s: 200",
    "qc_min_valid_fraction: 0.8",
    "qc_max_longest_gap_s: 5.0",
    "qc_max_interpolated_fraction: 0.2",
    "movement_cutoff_cm_s: 5",
    "integration_period_frames: 5",
    "contact_distance_cm: 6",
    "body_exclusion_distance_cm: 1",
    "proximity_range_cm: [6, 10]",
    "require_orientation: false",
    "min_bout_s: 0.0",
    "max_gap_s: 0.0",
    "write_manifest: true",
    on_cup_lines,
    extra_lines
  ), config_path)

  list(root = root, config = config_path,
       output = file.path(root, "output", "HAB"))
}

test_that("the SocP batch script runs end to end on synthetic data", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  fixture <- make_socp_fixture(temp_run_dir("socp"))
  output <- run_assay_script("DLCA_SocP v.1.1.0.R", fixture$config)
  expect_script_succeeded(output, "SocP script")

  combined <- file.path(fixture$output, "combined_output.csv")
  expect_true(file.exists(combined))
  result <- utils::read.csv(combined, stringsAsFactors = FALSE)

  expect_equal(nrow(result), 1L)
  expect_equal(result$phase, "HAB")
  expect_equal(result$ID, "mouse-2")
  # 90 frames on the left stimulus at 30 fps.
  expect_equal(result$contactLeft, 3, tolerance = 1e-6)
  expect_equal(result$contactRight, 0)
  # Metadata "R" maps the left chamber to novel in this repository.
  expect_equal(result$contactNovel, result$contactLeft)
  expect_equal(result$coordinateUnit, "cm")
  expect_equal(result$totalTime, 10)
  expect_equal(result$validTime, 10)
  expect_true(result$qcPass)

  expect_true(file.exists(file.path(fixture$output, "tracking_qc.csv")))
  manifest <- yaml::yaml.load_file(file.path(fixture$output, "run_manifest.yaml"))
  expect_equal(manifest$extra$assay, "SocP")
  expect_equal(manifest$extra$phase, "HAB")
  expect_equal(manifest$config$contact_distance_cm, 6)
})

test_that("the SocP batch script removes time on top of a cup when configured", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  fixture <- make_socp_fixture(temp_run_dir("socp_on_cup"), climb_frames = 60)
  output <- run_assay_script("DLCA_SocP v.1.1.0.R", fixture$config)
  expect_script_succeeded(output, "SocP script with on-cup exclusion")

  result <- utils::read.csv(file.path(fixture$output, "combined_output.csv"), stringsAsFactors = FALSE)
  # 90 frames of floor contact count; the 60 frames on top of the cup do not.
  expect_equal(result$contactLeft, 3, tolerance = 1e-6)
  expect_equal(result$contactRemovedLeft, 2, tolerance = 1e-6)
  expect_equal(result$onCupLeft, 2, tolerance = 1e-6)
  expect_equal(result$onCupRight, 0)
  expect_equal(result$validTime, 10)

  episodes <- utils::read.csv(file.path(fixture$output, "on_cup_episodes.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(episodes), 1L)
  expect_equal(episodes$object, "socl")
  expect_equal(episodes$decision, "on_object")
  expect_equal(episodes$Code, "B001")
  manifest <- yaml::yaml.load_file(file.path(fixture$output, "run_manifest.yaml"))
  expect_true(manifest$config$on_cup_exclusion)
})

test_that("the SocP batch script leaves out reviewed windows and counts hidden time as contact", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  root <- temp_run_dir("socp_hidden")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  decisions <- file.path(root, "decisions.csv")
  writeLines(c("file,target,start_s,end_s,decision,reason,reviewer",
               "B001_socp,all,9,10,leave_out,the last second,test"), decisions)
  hidden_lines <- c("hidden_contact: true", "hidden_min_s: 0.2", "hidden_max_distance_cm: 10",
                    "hidden_nose_off_body_cm: 10", "hidden_min_nose_lost_share: 0.5")

  # Without the rule: 10 + 10 frames of contact with the right cup.
  plain <- make_socp_fixture(file.path(root, "plain"), hidden_frames = 60)
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", plain$config), "SocP script")
  before <- utils::read.csv(file.path(plain$output, "combined_output.csv"), stringsAsFactors = FALSE)
  expect_equal(before$contactRight, 20 / 30, tolerance = 1e-6)
  expect_equal(before$totalTime, 10)
  expect_equal(before$leftOutS, 0)

  fixture <- make_socp_fixture(file.path(root, "rule"), hidden_frames = 60,
                               extra_lines = c(paste0("qc_decisions_file: ", decisions), hidden_lines))
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", fixture$config), "SocP script with the rules")
  result <- utils::read.csv(file.path(fixture$output, "combined_output.csv"), stringsAsFactors = FALSE)
  # the 60 hidden frames are contact with the right cup; the left cup is unchanged
  expect_equal(result$contactRight, 80 / 30, tolerance = 1e-6)
  expect_equal(result$hiddenRight, 2, tolerance = 1e-6)
  expect_equal(result$hiddenFamiliar, result$hiddenRight)     # "R": the right cup is familiar
  expect_equal(result$contactLeft, before$contactLeft)
  # the last second is gone from every measure
  expect_equal(result$totalTime, 9, tolerance = 1e-6)
  expect_equal(result$leftOutS, 1, tolerance = 1e-6)

  episodes <- utils::read.csv(file.path(fixture$output, "hidden_episodes.csv"), stringsAsFactors = FALSE)
  expect_equal(episodes$decision, "hidden")
  expect_equal(episodes$target, "socr")
  manifest <- yaml::yaml.load_file(file.path(fixture$output, "run_manifest.yaml"))
  expect_equal(manifest$extra$qc_decisions[[1]]$md5, unname(tools::md5sum(decisions)))

  # hidden_as unobserved: the same stretch is neither contact nor valid time
  unobserved <- make_socp_fixture(file.path(root, "unobserved"), hidden_frames = 60,
                                  extra_lines = c(paste0("qc_decisions_file: ", decisions), hidden_lines,
                                                  "hidden_as: unobserved"))
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", unobserved$config), "SocP script, hidden as unobserved")
  u <- utils::read.csv(file.path(unobserved$output, "combined_output.csv"), stringsAsFactors = FALSE)
  expect_equal(u$contactRight, 20 / 30, tolerance = 1e-6)
  expect_equal(u$hiddenRight, 2, tolerance = 1e-6)
  expect_equal(u$validTime, result$validTime - 2, tolerance = 1e-6)
  expect_equal(u$contactLeft, result$contactLeft)
  expect_equal(yaml::yaml.load_file(file.path(unobserved$output, "run_manifest.yaml"))$config$hidden_as, "unobserved")
  # and it needs the hidden rule
  bare <- make_socp_fixture(file.path(root, "bare"), extra_lines = "hidden_as: unobserved")
  output <- run_assay_script("DLCA_SocP v.1.1.0.R", bare$config)
  expect_true(!is.null(attr(output, "status")) && attr(output, "status") != 0)
  expect_true(any(grepl("hidden_as needs hidden_contact", output)))
})

test_that("the SocP batch script applies the reviewed class of hidden frames", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  root <- temp_run_dir("socp_hidden_review")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  hidden_lines <- c("hidden_contact: true", "hidden_min_s: 0.2", "hidden_max_distance_cm: 10",
                    "hidden_nose_off_body_cm: 10", "hidden_min_nose_lost_share: 0.5")
  # The fixture's hidden stretch at the right cup is rows 151-210, clip frames 150-209.
  review <- file.path(root, "review.csv")
  writeLines(c("file,target,first_frame,last_frame,class,evidence,reviewer",
               "B001_socp,socr,150,179,no_contact,head out of the zone,test",
               "B001_socp,socr,180,209,unobserved,tail only,test"), review)
  rule <- make_socp_fixture(file.path(root, "rule"), hidden_frames = 60, extra_lines = hidden_lines)
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", rule$config), "SocP script, hidden rule")
  r <- utils::read.csv(file.path(rule$output, "combined_output.csv"), stringsAsFactors = FALSE)
  fixture <- make_socp_fixture(file.path(root, "reviewed"), hidden_frames = 60,
                               extra_lines = c(hidden_lines, paste0("hidden_review_file: ", review)))
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", fixture$config), "SocP script, reviewed")
  x <- utils::read.csv(file.path(fixture$output, "combined_output.csv"), stringsAsFactors = FALSE)
  expect_equal(r$contactRight - x$contactRight, 2, tolerance = 1e-6)
  expect_equal(x$hiddenRight, 0)
  expect_equal(c(x$hiddenNoContact, x$hiddenUnobserved), c(1, 1), tolerance = 1e-6)
  expect_equal(x$validTime, r$validTime - 1, tolerance = 1e-6)
  expect_equal(x$contactLeft, r$contactLeft)
  episodes <- utils::read.csv(file.path(fixture$output, "hidden_episodes.csv"), stringsAsFactors = FALSE)
  expect_equal(episodes$decision, "hidden")
  expect_equal(unname(unlist(episodes[, c("review_contact_s", "review_no_contact_s", "review_unobserved_s",
                                          "unreviewed_s")])), c(0, 1, 1, 0), tolerance = 1e-6)
  manifest <- yaml::yaml.load_file(file.path(fixture$output, "run_manifest.yaml"))
  expect_equal(manifest$extra$hidden_review[[1]]$md5, unname(tools::md5sum(review)))

  # it needs the hidden rule
  bare <- make_socp_fixture(file.path(root, "bare"), extra_lines = paste0("hidden_review_file: ", review))
  output <- run_assay_script("DLCA_SocP v.1.1.0.R", bare$config)
  expect_true(!is.null(attr(output, "status")) && attr(output, "status") != 0)
  expect_true(any(grepl("hidden_review_file needs hidden_contact", output)))
})

test_that("the SocP batch script leaves out reviewed hands and clears a reviewed not-hidden stretch", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  root <- temp_run_dir("socp_hands")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  # A hand in view at the start until clip frame 29, and a hand that ends the phase at clip frame 270.
  hands <- file.path(root, "hands.csv")
  writeLines(c("file,kind,first_frame,last_frame,evidence,reviewer",
               "B001_socp,start,0,29,hand closing a cup,test",
               "B001_socp,end,270,,lid lifted,test"), hands)
  # The fixture's hidden stretch at the right cup is clip frames 151-210 (5.0 to 7.0 s), reviewed as not hidden;
  # counted on the rows after the start hand it would run 4.0-6.0 s, beyond the 0.5 s an edge may move.
  decisions <- file.path(root, "decisions.csv")
  writeLines(c("file,target,start_s,end_s,decision,reason,reviewer",
               "B001_socp,socr,5,7,not_hidden,mouse in view,test"), decisions)
  hidden_lines <- c("hidden_contact: true", "hidden_min_s: 0.2", "hidden_max_distance_cm: 10",
                    "hidden_nose_off_body_cm: 10", "hidden_min_nose_lost_share: 0.5")
  fixture <- make_socp_fixture(file.path(root, "rule"), hidden_frames = 60,
                               extra_lines = c(paste0("qc_decisions_file: ", decisions),
                                               paste0("hand_events_file: ", hands), hidden_lines))
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", fixture$config), "SocP script with hands")
  result <- utils::read.csv(file.path(fixture$output, "combined_output.csv"), stringsAsFactors = FALSE)
  # 30 frames at the start and 30 at the end are gone; contact with the left cup was frames 1-90
  expect_equal(result$leftOutS, 2, tolerance = 1e-6)
  expect_equal(result$totalTime, 8, tolerance = 1e-6)
  expect_equal(result$contactLeft, 2, tolerance = 1e-6)
  # the reviewed stretch is not contact: the right cup keeps its 20 seen frames
  expect_equal(result$hiddenRight, 0)
  expect_equal(result$contactRight, 20 / 30, tolerance = 1e-6)
  episodes <- utils::read.csv(file.path(fixture$output, "hidden_episodes.csv"), stringsAsFactors = FALSE)
  expect_equal(episodes$decision, "not_hidden")
  expect_equal(c(episodes$start_frame, episodes$start_s), c(151, 5))   # still on the clip's clock
  manifest <- yaml::yaml.load_file(file.path(fixture$output, "run_manifest.yaml"))
  expect_equal(manifest$extra$hand_events[[1]]$md5, unname(tools::md5sum(hands)))

  # a configured hand events file must exist
  writeLines(sub("^hand_events_file: .*$", paste0("hand_events_file: ", file.path(root, "absent.csv")),
                 readLines(fixture$config)), fixture$config)
  output <- run_assay_script("DLCA_SocP v.1.1.0.R", fixture$config)
  expect_true(!is.null(attr(output, "status")) && attr(output, "status") != 0)
  expect_true(any(grepl("hand_events_file not found", output)))
})

test_that("after a leave-out at the start, the SocP script keeps decisions and episode tables on the clip's clock", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  root <- temp_run_dir("socp_clip_clock")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  read_out <- function(fixture, name) utils::read.csv(file.path(fixture$output, name), stringsAsFactors = FALSE)
  header <- "file,target,start_s,end_s,decision,reason,reviewer"
  start_out <- "B001_socp,all,0,0.07,leave_out,the first three frames,test"

  # On-cup: the episode of a run without leave-out is the reference.
  ref <- make_socp_fixture(file.path(root, "cup_ref"), climb_frames = 60)
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", ref$config), "SocP script")
  ref_ep <- read_out(ref, "on_cup_episodes.csv")
  expect_equal(nrow(ref_ep), 1L)
  # A reviewed window that starts 0.45 s into the episode (clip time) still clears it; counted on the cut rows,
  # the episode would start 0.55 s before the window, beyond the 0.5 s an edge may move, and the run would stop.
  decisions <- file.path(root, "cup_decisions.csv")
  writeLines(c(header, start_out, sprintf("B001_socp,socl,%.4f,%.4f,not_on_object,reviewed,test",
                                          ref_ep$start_s + 0.45, ref_ep$start_s + ref_ep$duration_s)), decisions)
  cut <- make_socp_fixture(file.path(root, "cup_cut"), climb_frames = 60,
                           extra_lines = paste0("qc_decisions_file: ", decisions))
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", cut$config), "SocP script, start leave-out")
  cut_ep <- read_out(cut, "on_cup_episodes.csv")
  expect_equal(cut_ep[, c("start_frame", "end_frame", "start_s")], ref_ep[, c("start_frame", "end_frame", "start_s")])
  result <- read_out(cut, "combined_output.csv")
  expect_equal(result$leftOutS, 0.1, tolerance = 1e-6)
  expect_equal(result$onCupLeft, 0)

  # Hidden: the same stretch with and without the first three frames.
  hidden_lines <- c("hidden_contact: true", "hidden_min_s: 0.2", "hidden_max_distance_cm: 10",
                    "hidden_nose_off_body_cm: 10", "hidden_min_nose_lost_share: 0.5")
  ref <- make_socp_fixture(file.path(root, "hidden_ref"), hidden_frames = 60, extra_lines = hidden_lines)
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", ref$config), "SocP script, hidden")
  decisions <- file.path(root, "hidden_decisions.csv")
  writeLines(c(header, start_out), decisions)
  cut <- make_socp_fixture(file.path(root, "hidden_cut"), hidden_frames = 60,
                           extra_lines = c(paste0("qc_decisions_file: ", decisions), hidden_lines))
  expect_script_succeeded(run_assay_script("DLCA_SocP v.1.1.0.R", cut$config), "SocP script, hidden, start leave-out")
  cols <- c("target", "start_frame", "end_frame", "start_s", "end_s", "duration_s", "decision")
  expect_equal(read_out(cut, "hidden_episodes.csv")[, cols], read_out(ref, "hidden_episodes.csv")[, cols])
})

test_that("the SocP batch script repairs a misplaced nose when configured", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  fixture <- make_socp_fixture(temp_run_dir("socp_repair"))
  # Replace the animal by a full skeleton walking a circle (10 px per cm); its nose sits on a far point for 30 frames.
  n <- 300; t <- seq_len(n) / n * 2 * pi; heading <- t + pi / 2
  shape <- rbind(nose = c(4, 0), leftEar = c(2.5, 0.8), rightEar = c(2.5, -0.8), spine1 = c(1.5, 0), leftArm = c(1.5, 1),
                 rightArm = c(1.5, -1), bodycentre = c(0, 0), leftFlank = c(0, 1.4), rightFlank = c(0, -1.4), spine2 = c(-1.5, 0), tailbase = c(-2.7, 0))
  set.seed(3)
  points <- list(tl = cbind(x = rep(0, n), y = rep(240, n)), tr = cbind(x = rep(440, n), y = rep(240, n)),
                 br = cbind(x = rep(440, n), y = rep(0, n)), bl = cbind(x = rep(0, n), y = rep(0, n)),
                 socl = cbind(x = rep(50, n), y = rep(120, n)), socr = cbind(x = rep(390, n), y = rep(120, n)))
  for (p in rownames(shape)) points[[p]] <- 10 * cbind(x = 22 + 8 * cos(t) + shape[p, 1] * cos(heading) - shape[p, 2] * sin(heading) + stats::rnorm(n, 0, 0.03),
                                                       y = 12 + 8 * sin(t) + shape[p, 1] * sin(heading) + shape[p, 2] * cos(heading) + stats::rnorm(n, 0, 0.03))
  points$nose[101:130, ] <- cbind(x = 20, y = 220)
  write_tracking_fixture(file.path(fixture$root, "formatted", "HAB", "B001_socp.csv"), points, n)
  writeLines(c(readLines(fixture$config), "landmark_repair: true",
               "landmark_repair_landmarks: [nose, leftEar, rightEar, spine1, leftArm, rightArm, bodycentre, leftFlank, rightFlank, spine2, tailbase]",
               "landmark_repair_targets: [nose, bodycentre]", "landmark_repair_margin_cm: 3", "landmark_repair_margin_fraction: 0.6",
               "landmark_repair_min_share: 0.5", "landmark_repair_min_landmarks: 4", "landmark_repair_max_anchor_error_cm: 2",
               "landmark_repair_outside_cm: 15"), fixture$config)

  output <- run_assay_script("DLCA_SocP v.1.1.0.R", fixture$config)
  expect_script_succeeded(output, "SocP script with landmark repair")
  result <- utils::read.csv(file.path(fixture$output, "combined_output.csv"), stringsAsFactors = FALSE)
  expect_equal(result$noseReconstructedS, 1, tolerance = 1e-6)
  expect_equal(result$landmarksFlaggedS, 1, tolerance = 1e-6)
  expect_equal(result$skeletonCollapseS, 0)
  expect_lt(result$noseAnchorErrorP95, 0.5)
  expect_equal(result$validTime, 10)
  for (f in c("landmark_repairs.csv", "landmark_repair_episodes.csv", "landmark_repair_models.csv")) expect_true(file.exists(file.path(fixture$output, f)))
  ep <- utils::read.csv(file.path(fixture$output, "landmark_repair_episodes.csv"), stringsAsFactors = FALSE)
  expect_equal(ep$target, "nose")
  expect_equal(ep$reconstructed_share, 1)
  qc <- utils::read.csv(file.path(fixture$output, "tracking_qc.csv"), stringsAsFactors = FALSE)
  expect_equal(qc$reconstructed_time_s, 1, tolerance = 1e-6)
})
