testthat::test_that("SocInt rejects incompatible threshold units", {
  socint <- new.env(parent = globalenv())
  old_wd <- getwd()
  old_skip <- Sys.getenv("SLEAP_ANALYZER_SKIP_BATCH", unset = NA_character_)
  on.exit({
    setwd(old_wd)
    if (is.na(old_skip)) Sys.unsetenv("SLEAP_ANALYZER_SKIP_BATCH") else
      Sys.setenv(SLEAP_ANALYZER_SKIP_BATCH = old_skip)
  }, add = TRUE)
  setwd(file.path(repo_root, "02_SLEAPanalzyer"))
  Sys.setenv(SLEAP_ANALYZER_SKIP_BATCH = "true")
  sys.source("DLCA_SocInt v.0.0.2.r", envir = socint)

  invalid <- socint$config
  invalid$use_arena_calibration <- TRUE
  testthat::expect_error(socint$validate_socint_units(invalid), "does not match")
  invalid$threshold_unit <- "cm"
  testthat::expect_error(socint$validate_socint_units(invalid), "not been scientifically confirmed")
  invalid$scientific_cm_thresholds_confirmed <- TRUE
  testthat::expect_silent(socint$validate_socint_units(invalid))
  testthat::expect_silent(socint$validate_socint_units(socint$config))
})

testthat::test_that("SocInt event helpers match required timing behavior", {
  socint <- new.env(parent = globalenv())
  old_wd <- getwd()
  old_skip <- Sys.getenv("SLEAP_ANALYZER_SKIP_BATCH", unset = NA_character_)
  on.exit({
    setwd(old_wd)
    if (is.na(old_skip)) Sys.unsetenv("SLEAP_ANALYZER_SKIP_BATCH") else
      Sys.setenv(SLEAP_ANALYZER_SKIP_BATCH = old_skip)
  }, add = TRUE)
  setwd(file.path(repo_root, "02_SLEAPanalzyer"))
  Sys.setenv(SLEAP_ANALYZER_SKIP_BATCH = "true")
  sys.source("DLCA_SocInt v.0.0.2.r", envir = socint)

  event <- c(rep(FALSE, 60), rep(TRUE, 30), rep(FALSE, 60))
  result <- socint$summarise_event(event, fps = 30)
  testthat::expect_equal(result$duration_s, 1)
  testthat::expect_equal(result$bout_n, 1)
  testthat::expect_equal(result$latency_s, 2)
  testthat::expect_true(is.na(socint$summarise_event(rep(FALSE, 10), 30)$latency_s))

  short_gap <- socint$interpolate_with_qc(c(NA, 1, 2, NA), maxgap = 1)$x
  testthat::expect_equal(short_gap, c(1, 1, 2, 2))
  long_gap <- socint$interpolate_with_qc(c(NA, NA, 1, 2), maxgap = 1)$x
  testthat::expect_true(all(is.na(long_gap[1:2])))
})

testthat::test_that("SocInt arena calibration reports anisotropy and rejects degeneration", {
  socint <- new.env(parent = globalenv())
  old_wd <- getwd()
  old_skip <- Sys.getenv("SLEAP_ANALYZER_SKIP_BATCH", unset = NA_character_)
  on.exit({
    setwd(old_wd)
    if (is.na(old_skip)) Sys.unsetenv("SLEAP_ANALYZER_SKIP_BATCH") else
      Sys.setenv(SLEAP_ANALYZER_SKIP_BATCH = old_skip)
  }, add = TRUE)
  setwd(file.path(repo_root, "02_SLEAPanalzyer"))
  Sys.setenv(SLEAP_ANALYZER_SKIP_BATCH = "true")
  sys.source("DLCA_SocInt v.0.0.2.r", envir = socint)

  geom <- tibble::tibble(
    corner = c("tl", "tr", "br", "bl"),
    x = c(0, 100, 100, 0), y = c(100, 100, 0, 0)
  )
  config <- socint$config
  config$arena_width_cm <- 50
  config$arena_height_cm <- 50
  calibration <- socint$estimate_arena_calibration(geom, config)
  testthat::expect_equal(calibration$cm_per_px_x, 0.5)
  testthat::expect_equal(calibration$cm_per_px_y, 0.5)
  testthat::expect_equal(calibration$calibration_anisotropy_percent, 0)

  config$arena_height_cm <- 25
  anisotropic <- socint$estimate_arena_calibration(geom, config)
  testthat::expect_gt(anisotropic$calibration_anisotropy_percent, 0)

  degenerate <- geom
  degenerate$x <- 0
  degenerate$y <- 0
  testthat::expect_error(socint$estimate_arena_calibration(degenerate, config), "degenerate")
})
