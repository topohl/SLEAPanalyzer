# End-to-end run of the SocP batch script against synthetic data.

make_socp_fixture <- function(root, n = 300, fps = 30, contact_frames = 90,
                              climb_frames = 0) {
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
    on_cup_lines
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
