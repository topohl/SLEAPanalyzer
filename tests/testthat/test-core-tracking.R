testthat::test_that("DLC/SLEAP CSV parsing handles missing coordinates and frame origins", {
  csv <- tempfile(fileext = ".csv")
  writeLines(c(
    "scorer,model,model,model,model,model,model",
    "bodyparts,nose,nose,nose,bodycentre,bodycentre,bodycentre",
    "coords,x,y,likelihood,x,y,likelihood",
    "100,1,2,0.9,10,20,1",
    "101,,4,0.8,12,22,1",
    "102,5,6,0.7,14,24,1"
  ), csv)

  tracking <- ReadDLCDataFromCSV(csv, fps = 10)
  testthat::expect_equal(tracking$frames, 100:102)
  testthat::expect_equal(tracking$seconds, c(0, 0.1, 0.2))
  testthat::expect_equal(tracking$median.data["nose", "x"], 3)
  testthat::expect_equal(tracking$filename, basename(csv))
})

testthat::test_that("trimming removes exact counts and keeps every component aligned", {
  tracking <- make_tracking(
    frames = 100:109, fps = 10,
    points = list(
      bodycentre = cbind(x = 1:10, y = 11:20),
      nose = cbind(x = 21:30, y = 31:40)
    )
  )
  tracking$labels$manual <- letters[1:10]
  tracking$features <- data.frame(feature = 101:110)

  start_cut <- CutTrackingData(tracking, start = 2)
  testthat::expect_equal(start_cut$frames, 102:109)
  testthat::expect_equal(start_cut$seconds, tracking$seconds[3:10])
  testthat::expect_equal(start_cut$labels$manual, letters[3:10])
  testthat::expect_equal(start_cut$features$feature, 103:110)
  testthat::expect_equal(start_cut$data$nose$frame, 102:109)

  end_cut <- CutTrackingData(tracking, end = 2)
  testthat::expect_equal(end_cut$frames, 100:107)
  testthat::expect_equal(nrow(end_cut$data$bodycentre), 8)

  selected <- CutTrackingData(tracking, remove.frames = c(102, 105), keep.frames = c(101, 103, 109))
  testthat::expect_equal(selected$frames, c(101, 103, 109))
  testthat::expect_equal(selected$data$bodycentre$x, c(2, 4, 10))

  binned <- AddBinData(
    tracking,
    bindat = data.frame(bin = "window", from = 0.2, to = 0.5),
    unit = "second"
  )
  testthat::expect_equal(binned$bins$from, 102)
  testthat::expect_equal(binned$bins$to, 105)
})

testthat::test_that("movement is zero when stationary and correct at constant velocity", {
  stationary <- make_tracking(
    frames = 0:30, fps = 10,
    points = list(bodycentre = cbind(x = rep(5, 31), y = rep(-2, 31)))
  )
  stationary <- CalculateMovement(stationary, movement_cutoff = 1, integration_period = 0)
  testthat::expect_equal(stationary$data$bodycentre$speed, rep(0, 31))
  testthat::expect_equal(sum(stationary$data$bodycentre$speed), 0)

  moving <- make_tracking(
    frames = 0:30, fps = 10,
    points = list(bodycentre = cbind(x = 2 * (0:30), y = rep(0, 31)))
  )
  moving <- CalculateMovement(moving, movement_cutoff = 1, integration_period = 0)
  testthat::expect_equal(sum(moving$data$bodycentre$speed), 60)
  testthat::expect_equal(moving$data$bodycentre$speed[-1] * moving$fps, rep(20, 30))
})

testthat::test_that("area calibration is correct, invariant, and rejects invalid geometry", {
  points <- square_points()
  points$bodycentre <- cbind(x = rep(20, 10), y = rep(40, 10))
  calibrated <- CalibrateTrackingData(
    make_tracking(points = points), method = "area", in.metric = 50 * 50,
    points = c("tl", "tr", "br", "bl")
  )
  testthat::expect_equal(calibrated$px.to.cm, 0.5)
  testthat::expect_equal(calibrated$distance.units, "cm")
  testthat::expect_equal(calibrated$data$bodycentre$x, rep(10, 10))

  translated_points <- square_points(offset = c(700, -300))
  translated_points$bodycentre <- cbind(x = rep(720, 10), y = rep(-260, 10))
  translated <- CalibrateTrackingData(
    make_tracking(points = translated_points), method = "area", in.metric = 2500,
    points = c("tl", "tr", "br", "bl")
  )
  testthat::expect_equal(
    translated$data$bodycentre$x - translated$data$tl$x,
    calibrated$data$bodycentre$x - calibrated$data$tl$x
  )

  scaled_points <- square_points(scale = 2)
  scaled_points$bodycentre <- cbind(x = rep(40, 10), y = rep(80, 10))
  scaled <- CalibrateTrackingData(
    make_tracking(points = scaled_points), method = "area", in.metric = 2500,
    points = c("tl", "tr", "br", "bl")
  )
  testthat::expect_equal(scaled$data$bodycentre$x, calibrated$data$bodycentre$x)
  testthat::expect_error(CalibrateTrackingData(calibrated, "ratio", ratio = 1), "already calibrated")

  degenerate <- square_points()
  degenerate$tr <- degenerate$tl
  degenerate$br <- degenerate$bl
  testthat::expect_error(
    CalibrateTrackingData(make_tracking(points = degenerate), "area", 2500, c("tl", "tr", "br", "bl")),
    "Degenerate"
  )
  testthat::expect_error(CalibrateTrackingData(make_tracking(points = points), "ratio", ratio = 0), "positive")
})

testthat::test_that("OFT zones require four valid corners and classify center/periphery", {
  testthat::skip_if_not_installed("sp")
  points <- square_points(n = 2)
  points$bodycentre <- cbind(x = c(50, 5), y = c(50, 50))
  tracking <- AddOFTZones(make_tracking(frames = 0:1, points = points))

  testthat::expect_true(IsInZone(tracking, "bodycentre", "center")[1])
  testthat::expect_true(IsInZone(tracking, "bodycentre", "periphery", invert = TRUE)[2])
  testthat::expect_error(AddOFTZones(make_tracking(frames = 0:1, points = points), c("tl", "tr", "br")), "Exactly four")
  testthat::expect_error(AddOFTZones(make_tracking(frames = 0:1, points = points), c("tl", "tr", "br", "missing")), "not found")
})
