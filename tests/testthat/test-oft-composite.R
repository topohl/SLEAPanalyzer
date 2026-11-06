# The OFT composite score is cohort-relative, not a subject-level measurement.

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tidyr")
suppressPackageStartupMessages(library(dplyr))

oft_functions <- new.env(parent = environment())
local({
  script <- readLines(
    file.path(repo_root, "02_SLEAPanalzyer", "DLCA_OFT v1.2.0.R"), warn = FALSE
  )
  # Evaluate only the helper definitions, not the batch loop.
  start <- grep("^zscore <- function", script)[1]
  stop_at <- grep("^add_center_exploration_score <- function", script)[1]
  end <- grep("^make_theme_figure <- function", script)[1] - 1L
  eval(parse(text = paste(script[c(start:(start + 6), stop_at:end)], collapse = "\n")),
       envir = oft_functions)
})

cohort <- function(center_time_percent, center_entries, center_latency_s,
                   mean_wall_distance_cm) {
  data.frame(
    center_time_percent = center_time_percent,
    center_entries = center_entries,
    center_latency_s = center_latency_s,
    mean_wall_distance_cm = mean_wall_distance_cm
  )
}

test_that("the composite is standardised within the analysed cohort", {
  scored <- oft_functions$add_center_exploration_score(
    cohort(c(10, 20, 30), c(2, 4, 6), c(60, 30, 10), c(20, 15, 10))
  )
  # z-scores of a cohort sum to zero by construction.
  expect_equal(sum(scored$z_center_time_percent), 0, tolerance = 1e-12)
  expect_equal(scored$oft_center_exploration_cohort_n, rep(3L, 3))
})

test_that("the same animal scores differently in a different cohort", {
  # One animal, identical raw measurements, analysed with two different
  # cohorts. If the score were a subject-level measurement it would not move.
  target <- c(20, 4, 30, 15)

  low_cohort <- oft_functions$add_center_exploration_score(
    cohort(c(target[1], 5, 6), c(target[2], 1, 1), c(target[3], 90, 95), c(target[4], 22, 23))
  )
  high_cohort <- oft_functions$add_center_exploration_score(
    cohort(c(target[1], 50, 60), c(target[2], 12, 14), c(target[3], 2, 3), c(target[4], 4, 5))
  )

  low <- low_cohort$oft_center_exploration_score_cohort_z_experimental[1]
  high <- high_cohort$oft_center_exploration_score_cohort_z_experimental[1]
  expect_false(isTRUE(all.equal(low, high)))
  # Against a low-exploring cohort the animal looks high, and vice versa.
  expect_gt(low, high)
})

test_that("the raw components are retained alongside the composite", {
  scored <- oft_functions$add_center_exploration_score(
    cohort(c(10, 20, 30), c(2, 4, 6), c(60, 30, 10), c(20, 15, 10))
  )
  # Modelling the raw components is preferable to modelling a composite whose
  # scale depends on the cohort.
  expect_true(all(c(
    "center_time_percent", "center_entries", "center_latency_s",
    "mean_wall_distance_cm"
  ) %in% names(scored)))
})

test_that("the composite column name carries its caveats", {
  scored <- oft_functions$add_center_exploration_score(
    cohort(c(10, 20, 30), c(2, 4, 6), c(60, 30, 10), c(20, 15, 10))
  )
  expect_true(
    "oft_center_exploration_score_cohort_z_experimental" %in% names(scored)
  )
  # The bare name implied an intrinsic subject-level property.
  expect_false("oft_center_exploration_score" %in% names(scored))
})

test_that("a cohort with no variation yields NA rather than a fabricated score", {
  scored <- oft_functions$add_center_exploration_score(
    cohort(c(20, 20, 20), c(4, 4, 4), c(30, 30, 30), c(15, 15, 15))
  )
  expect_true(all(is.na(scored$oft_center_exploration_score_cohort_z_experimental)))
})

test_that("a partial last bin does not enter the habituation slope", {
  bins_env <- new.env(parent = environment())
  script <- readLines(
    file.path(repo_root, "02_SLEAPanalzyer", "DLCA_OFT v1.2.0.R"), warn = FALSE
  )
  start <- grep("^make_binned_table <- function", script)[1]
  end <- grep("^make_epoch_metrics <- function", script)[1] - 1L
  eval(parse(text = paste(script[start:end], collapse = "\n")), envir = bins_env)
  testthat::skip_if(!exists("make_habituation_slopes", envir = bins_env), "slope helpers not found")
  testthat::skip_if_not_installed("purrr")
  config <- list(bin_size_s = 60, fps = 30)
  # 10 full minutes with distance falling by 10 cm per minute, then 5 frames of an 11th minute
  n <- 10 * 1800 + 5
  time_s <- (seq_len(n) - 1) / 30
  minute <- floor(time_s / 60)
  frame_tbl <- data.frame(file = "f", Batch = "B1", Code = "AAAA", time_s = time_s,
                          in_center = FALSE, in_periphery = TRUE, in_corners = FALSE,
                          is_moving = TRUE, is_immobile = FALSE, wall_distance_cm = 5,
                          speed_cm_s = (200 - 10 * minute) / 60)
  bins_env$entry_count <- function(x) 0L
  binned <- bins_env$make_binned_table(frame_tbl, config)
  expect_equal(nrow(binned), 11L)
  expect_equal(binned$complete_bin, c(rep(TRUE, 10), FALSE))
  # the synthetic minutes fit exactly, which lm() warns about
  slopes <- suppressWarnings(bins_env$make_habituation_slopes(binned))
  expect_equal(slopes$slope_per_min_distance_cm, -10, tolerance = 1e-6)
})

test_that("habituation p-values are labelled as within-animal and nominal", {
  script <- readLines(
    file.path(repo_root, "02_SLEAPanalzyer", "DLCA_OFT v1.2.0.R"), warn = FALSE
  )
  # The regression treats one animal's consecutive time bins as independent
  # observations, which they are not; the column name must say so.
  expect_true(any(grepl("p_value_within_animal_nominal", script, fixed = TRUE)))
  expect_false(any(grepl("^    p_value = sm\\$coefficients", script)))
})
