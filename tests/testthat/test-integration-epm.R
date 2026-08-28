# End-to-end run of the EPM batch script using the shipped zone definitions.
#
# Maze in pixels, 10 px per cm over a 60 x 60 cm extent. Arms are 10 cm wide.
#
#                 open.top          y 350..600
#   closed.left     center     closed.right    y 250..350
#                 open.bottom       y 0..250

epm_landmarks <- function(n) {
  at <- function(x, y) cbind(x = rep(x, n), y = rep(y, n))
  list(
    tl = at(250, 600), tr = at(350, 600),
    ctr = at(350, 350), ctl = at(250, 350),
    cbl = at(250, 250), cbr = at(350, 250),
    br = at(350, 0), bl = at(250, 0),
    lt = at(0, 350), lb = at(0, 250),
    rt = at(600, 350), rb = at(600, 250)
  )
}

make_epm_fixture <- function(root, fps = 10) {
  # 40 frames: 10 in the left closed arm, 10 in the centre, 10 in the top open
  # arm, 10 back in the centre.
  segments <- list(
    c(x = 120, y = 300),  # closed left
    c(x = 300, y = 300),  # centre
    c(x = 300, y = 500),  # open top
    c(x = 300, y = 300)   # centre
  )
  body_x <- unlist(lapply(segments, function(s) rep(s[["x"]], 10)))
  body_y <- unlist(lapply(segments, function(s) rep(s[["y"]], 10)))
  n <- length(body_x)

  points <- epm_landmarks(n)
  points$bodycentre <- cbind(x = body_x, y = body_y)
  points$headcentre <- cbind(x = body_x, y = body_y)
  points$neck <- cbind(x = body_x, y = body_y)

  write_tracking_fixture(file.path(root, "formatted", "E001_epm.csv"), points, n)

  config_path <- file.path(root, "epm.yaml")
  writeLines(c(
    paste0("input_dir: ", file.path(root, "formatted")),
    paste0("output_dir: ", file.path(root, "output")),
    paste0("zone_file: ", file.path(repo_root, "02_SLEAPanalzyer", "EPM_zoneinfo.csv")),
    paste0("fps: ", fps),
    "arena_width_cm: 60",
    "arena_height_cm: 60",
    "arena_corner_names: [tl, tr, br, bl]",
    "max_interpolation_gap_s: 0.2",
    "likelihood_cutoff: ~",
    "max_plausible_speed_cm_s: 500",
    "qc_min_valid_fraction: 0.8",
    "qc_max_longest_gap_s: 5.0",
    "qc_max_interpolated_fraction: 0.2",
    "movement_cutoff_cm_s: 5",
    "integration_period_frames: 0",
    "nose_dips: false",
    "min_bout_s: 0.0",
    "max_gap_s: 0.0",
    "write_manifest: true"
  ), config_path)

  list(root = root, config = config_path, output = file.path(root, "output"))
}

test_that("the EPM batch script runs end to end and assigns arms correctly", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("cowplot")

  fixture <- make_epm_fixture(temp_run_dir("epm"))
  output <- run_assay_script("DLCA_EPM v1.0.0.R", fixture$config)
  expect_script_succeeded(output, "EPM script")

  report_path <- file.path(fixture$output, "Report.csv")
  expect_true(file.exists(report_path))
  report <- utils::read.csv(report_path, stringsAsFactors = FALSE)
  expect_equal(nrow(report), 1L)

  # 10 frames per segment at 10 fps is 1 s each.
  expect_equal(report$bodycentre.closed.total.time, 1)
  expect_equal(report$bodycentre.closed.left.total.time, 1)
  expect_equal(report$bodycentre.closed.right.total.time, 0)
  expect_equal(report$bodycentre.open.total.time, 1)
  expect_equal(report$bodycentre.open.top.total.time, 1)
  expect_equal(report$bodycentre.center.total.time, 2)
  expect_equal(report$bodycentre.total.time, 4)
  expect_equal(report$bodycentre.valid.time, 4)

  # One observed entry into the open arm; the opening centre visit has no
  # observed onset.
  expect_equal(report$bodycentre.open.entries, 1)
  expect_equal(report$bodycentre.closed.entries, 0)

  expect_true(file.exists(file.path(fixture$output, "tracking_qc.csv")))
  manifest <- yaml::yaml.load_file(file.path(fixture$output, "run_manifest.yaml"))
  expect_equal(manifest$extra$assay, "EPM")
})

test_that("the shipped EPM zone definitions are simple, non-degenerate polygons", {
  zone_file <- file.path(repo_root, "02_SLEAPanalzyer", "EPM_zoneinfo.csv")
  zones <- utils::read.table(
    zone_file, sep = ";", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
  )
  landmarks <- epm_landmarks(1)
  for (zone in names(zones)) {
    names_in_zone <- zones[[zone]]
    names_in_zone <- names_in_zone[!is.na(names_in_zone) & nzchar(trimws(names_in_zone))]
    polygon <- do.call(rbind, lapply(names_in_zone, function(nm) {
      data.frame(x = landmarks[[nm]][1, "x"], y = landmarks[[nm]][1, "y"])
    }))
    expect_gt(polygon_area(polygon), 0)
    expect_false(is_self_intersecting(polygon), info = zone)
  }
})

test_that("AddZones refuses a zone whose corners collapse to zero area", {
  n <- 5
  points <- epm_landmarks(n)
  points$bodycentre <- cbind(x = rep(300, n), y = rep(300, n))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)

  # Corners listed diagonally rather than around the perimeter. sp would
  # report almost every interior frame as outside this shape, silently
  # under-counting occupancy, so building it must fail instead.
  crossed <- data.frame(bad = c("lt", "cbl", "ctl", "lb"), stringsAsFactors = FALSE)
  expect_error(AddZones(tracking, crossed), "degenerate")
})

test_that("AddZones refuses a self-intersecting zone with non-zero area", {
  n <- 5
  points <- epm_landmarks(n)
  points$bodycentre <- cbind(x = rep(300, n), y = rep(300, n))
  # An asymmetric bowtie: the shoelace area does not cancel to zero, so only
  # the explicit self-intersection test catches it.
  points$p1 <- cbind(x = rep(0, n), y = rep(0, n))
  points$p2 <- cbind(x = rep(100, n), y = rep(0, n))
  points$p3 <- cbind(x = rep(20, n), y = rep(100, n))
  points$p4 <- cbind(x = rep(80, n), y = rep(-30, n))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)

  bowtie <- data.frame(bad = c("p1", "p2", "p3", "p4"), stringsAsFactors = FALSE)
  expect_gt(polygon_area(tracking$median.data[c("p1", "p2", "p3", "p4"), c("x", "y")]), 0)
  expect_error(AddZones(tracking, bowtie), "self-intersecting")
})

test_that("MultiFileReport binds rows without an undeclared data.table dependency", {
  # MultiFileReport() called rbindlist() unqualified, so every EPM run failed
  # unless data.table happened to be attached.
  expect_false("rbindlist" %in% ls(envir = globalenv()))
  bound <- bind_report_rows(list(
    list(file = "a", open.time = 1, closed.time = 2),
    list(file = "b", open.time = 3)
  ))
  expect_equal(nrow(bound), 2L)
  # A missing entry is filled with NA rather than dropping the row.
  expect_true(is.na(bound$closed.time[2]))
})

test_that("AddZones reports zones referring to untracked landmarks", {
  n <- 5
  points <- epm_landmarks(n)
  points$bodycentre <- cbind(x = rep(300, n), y = rep(300, n))
  tracking <- make_tracking(frames = 0:(n - 1), fps = 10, points = points)

  bad <- data.frame(zone = c("lt", "ctl", "cbl", "nonexistent"), stringsAsFactors = FALSE)
  expect_error(AddZones(tracking, bad), "untracked landmark")
})
