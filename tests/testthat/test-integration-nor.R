# End-to-end run of the NOR batch script against synthetic data.
#
# This is the test that proves the whole pipeline executes: configuration
# loading, tracking import, bounded interpolation, calibration, QC, the shared
# event engine, output writing and the run manifest.

make_nor_fixture <- function(root, n = 300, fps = 30, contact_frames = 60) {
  formatted <- file.path(root, "formatted")
  dir.create(formatted, recursive = TRUE, showWarnings = FALSE)

  # Arena corners at 0..490 px for a 49 x 49 cm arena, so 10 px per cm.
  corner <- function(x, y) cbind(x = rep(x, n), y = rep(y, n))
  # Objects 15 cm from each side wall, mid height.
  objL <- corner(150, 245)
  objR <- corner(350, 245)

  # The nose sits on the left object for the first contact_frames frames, then
  # parks in the middle of the arena, far from both objects.
  nose_x <- c(rep(150, contact_frames), rep(250, n - contact_frames))
  nose_y <- c(rep(215, contact_frames), rep(100, n - contact_frames))
  # bodycentre behind the nose so the animal faces the object.
  body_x <- nose_x
  body_y <- nose_y - 60

  points <- list(
    tl = corner(0, 490), tr = corner(490, 490),
    br = corner(490, 0), bl = corner(0, 0),
    objL = objL, objR = objR,
    nose = cbind(x = nose_x, y = nose_y),
    bodycentre = cbind(x = body_x, y = body_y)
  )

  header_names <- c("scorer")
  header_coords <- c("bodyparts")
  body <- data.frame(frame = 0:(n - 1))
  for (nm in names(points)) {
    header_names <- c(header_names, rep(nm, 3))
    header_coords <- c(header_coords, "x", "y", "likelihood")
    body[[paste0(nm, "_x")]] <- points[[nm]][, "x"]
    body[[paste0(nm, "_y")]] <- points[[nm]][, "y"]
    body[[paste0(nm, "_l")]] <- rep(1, n)
  }
  lines <- c(
    paste(rep("scorer", length(header_names)), collapse = ","),
    paste(header_names, collapse = ","),
    paste(header_coords, collapse = ","),
    apply(body, 1, function(r) paste(r, collapse = ","))
  )
  writeLines(lines, file.path(formatted, "A001_nor.csv"))

  writeLines(c("Code\tNovelLoc", "A001\tR"), file.path(root, "novelLoc.txt"))
  writeLines(c("Code ID", "A001 mouse-1"), file.path(root, "animalIDCode.txt"))

  config_path <- file.path(root, "nor.yaml")
  writeLines(c(
    paste0("input_dir: ", file.path(root, "formatted")),
    paste0("output_dir: ", file.path(root, "output")),
    paste0("metadata_dir: ", root),
    paste0("animal_id_code_file: ", file.path(root, "animalIDCode.txt")),
    "novel_location_file: novelLoc.txt",
    paste0("fps: ", fps),
    "arena_width_cm: 49",
    "arena_height_cm: 49",
    "arena_corner_names: [tl, tr, br, bl]",
    "max_interpolation_gap_s: 0.2",
    "likelihood_cutoff: ~",
    "max_plausible_speed_cm_s: 200",
    "qc_min_valid_fraction: 0.8",
    "qc_max_longest_gap_s: 5.0",
    "qc_max_interpolated_fraction: 0.2",
    "movement_cutoff_cm_s: 5",
    "integration_period_frames: 5",
    "contact_geometry: radial",
    "contact_distance_cm: 4",
    "body_exclusion_distance_cm: 1",
    "object_box_width_cm: 9",
    "object_box_height_cm: 7",
    "contact_angle_deg: [70, 290]",
    "proximity_range_cm: [4, 8]",
    "proximity_angle_deg: [90, 270]",
    "min_bout_s: 0.0",
    "max_gap_s: 0.0",
    "report_rearing: false",
    "write_manifest: true"
  ), config_path)

  list(root = root, config = config_path, output = file.path(root, "output"))
}

run_nor_script <- function(fixture) {
  script <- file.path(repo_root, "02_SLEAPanalzyer", "DLCA_NOR v1.2.1.R")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  # system2(env = ) is ignored on Windows, so set the variable in this process
  # and let the child inherit it.
  previous <- Sys.getenv("SLEAP_ANALYZER_CONFIG", unset = NA_character_)
  Sys.setenv(SLEAP_ANALYZER_CONFIG = fixture$config)
  on.exit({
    if (is.na(previous)) Sys.unsetenv("SLEAP_ANALYZER_CONFIG") else
      Sys.setenv(SLEAP_ANALYZER_CONFIG = previous)
  }, add = TRUE)
  suppressWarnings(system2(
    rscript, c("--vanilla", shQuote(script)),
    stdout = TRUE, stderr = TRUE
  ))
}

test_that("the NOR batch script runs end to end on synthetic data", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  fixture <- make_nor_fixture(file.path(tempdir(), paste0("nor-", as.integer(runif(1, 1, 1e9)))))
  output <- run_nor_script(fixture)
  status <- attr(output, "status")

  if (!is.null(status) && status != 0) {
    fail(paste0("NOR script exited with status ", status, ":\n",
                paste(output, collapse = "\n")))
  }

  combined <- file.path(fixture$output, "combined_output.csv")
  expect_true(file.exists(combined))
  result <- utils::read.csv(combined, stringsAsFactors = FALSE)
  expect_equal(nrow(result), 1L)

  # 60 frames of contact with the left object at 30 fps.
  expect_equal(result$contactLeft, 2, tolerance = 1e-6)
  expect_equal(result$contactRight, 0)
  # Metadata "R" maps the left side to novel in this repository's convention.
  expect_equal(result$contactNov, result$contactLeft)
  expect_equal(result$contactFam, 0)
  expect_equal(result$novelLoc, "R")
  expect_equal(result$ID, "mouse-1")
  expect_equal(result$Code, "A001")
  expect_equal(result$contactGeometry, "radial")

  # Nothing was interpolated, so valid time is the whole recording.
  expect_equal(result$totalTime, 10)
  expect_equal(result$validTime, 10)
  expect_equal(result$validFraction, 1)
  expect_true(result$qcPass)
})

test_that("the run writes a QC table and a provenance manifest", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  fixture <- make_nor_fixture(file.path(tempdir(), paste0("nor-prov-", as.integer(runif(1, 1, 1e9)))))
  output <- run_nor_script(fixture)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    fail(paste0("NOR script exited with status ", status, ":\n",
                paste(output, collapse = "\n")))
  }

  qc_path <- file.path(fixture$output, "tracking_qc.csv")
  expect_true(file.exists(qc_path))
  qc <- utils::read.csv(qc_path, stringsAsFactors = FALSE)
  expect_equal(qc$valid_fraction, 1)
  expect_equal(qc$coordinate_unit, "cm")

  manifest_path <- file.path(fixture$output, "run_manifest.yaml")
  expect_true(file.exists(manifest_path))
  manifest <- yaml::yaml.load_file(manifest_path)

  expect_equal(manifest$extra$assay, "NOR")
  expect_equal(manifest$config$fps, 30)
  expect_equal(manifest$config$contact_distance_cm, 4)
  expect_match(manifest$run$timestamp_utc, "^\\d{4}-\\d{2}-\\d{2}T")
  expect_true(!is.null(manifest$environment$r_version))
  # The input file is identified by hash so a rerun is verifiable.
  expect_length(manifest$inputs, 1L)
  expect_equal(nchar(manifest$inputs[[1]]$md5), 32L)
})

test_that("an invalid configuration stops the run instead of guessing", {
  testthat::skip_if_not_installed("yaml")
  fixture <- make_nor_fixture(file.path(tempdir(), paste0("nor-bad-", as.integer(runif(1, 1, 1e9)))))
  # Corrupt one value.
  config_lines <- readLines(fixture$config)
  config_lines <- sub("^contact_geometry: radial$", "contact_geometry: teleport", config_lines)
  writeLines(config_lines, fixture$config)

  output <- run_nor_script(fixture)
  status <- attr(output, "status")
  expect_true(!is.null(status) && status != 0)
  expect_match(paste(output, collapse = "\n"), "contact_geometry")
})
