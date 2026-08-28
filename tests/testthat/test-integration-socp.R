# End-to-end run of the SocP batch script against synthetic data.

make_socp_fixture <- function(root, n = 300, fps = 30, contact_frames = 90) {
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
    "write_manifest: true"
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
