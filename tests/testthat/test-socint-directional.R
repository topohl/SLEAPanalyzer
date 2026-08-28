# End-to-end SocInt directional attribution on a synthetic dyad.

load_socint <- function() {
  socint <- new.env(parent = globalenv())
  old_wd <- getwd()
  old_skip <- Sys.getenv("SLEAP_ANALYZER_SKIP_BATCH", unset = NA_character_)
  on.exit({
    setwd(old_wd)
    if (is.na(old_skip)) Sys.unsetenv("SLEAP_ANALYZER_SKIP_BATCH") else
      Sys.setenv(SLEAP_ANALYZER_SKIP_BATCH = old_skip)
  }, add = TRUE, after = FALSE)
  setwd(file.path(repo_root, "02_SLEAPanalzyer"))
  Sys.setenv(SLEAP_ANALYZER_SKIP_BATCH = "true")
  sys.source("DLCA_SocInt v.0.0.2.r", envir = socint)
  socint
}

#' Build a two-animal Tracking object from body-centre trajectories.
#'
#' Each animal is a rigid body: nose one unit ahead of the body centre along
#' the direction of travel, tail base one unit behind, ears and sides offset
#' laterally. That keeps headings well defined without encoding any social
#' behavior beyond the supplied trajectories.
socint_tracking <- function(a_x, a_y, b_x, b_y, fps = 30) {
  n <- length(a_x)
  stopifnot(length(a_y) == n, length(b_x) == n, length(b_y) == n)

  animal <- function(cx, cy, suffix) {
    dx <- c(diff(cx), tail(diff(cx), 1))
    dy <- c(diff(cy), tail(diff(cy), 1))
    norm <- sqrt(dx^2 + dy^2)
    norm[norm == 0] <- 1
    ux <- dx / norm
    uy <- dy / norm
    # Left-hand normal of the heading.
    nx <- -uy
    ny <- ux
    at <- function(along, across) {
      data.frame(
        frame = seq_len(n) - 1L,
        x = cx + along * ux + across * nx,
        y = cy + along * uy + across * ny,
        likelihood = rep(1, n)
      )
    }
    out <- list(
      at(1, 0), at(1, 0.5), at(1, -0.5), at(0, 0),
      at(0, 0.5), at(0, -0.5), at(-1, 0), at(-2, 0)
    )
    names(out) <- paste0(
      c("nose", "leftEar", "rightEar", "bodycentre",
        "leftSide", "rightSide", "tailBase", "tailEnd"),
      suffix
    )
    out
  }

  data <- c(animal(a_x, a_y, "_1"), animal(b_x, b_y, "_2"))
  list(
    data = data,
    frames = seq_len(n) - 1L,
    fps = fps,
    seconds = (seq_len(n) - 1L) / fps,
    point.info = data.frame(PointName = names(data), PointType = "NotDefined"),
    distance.units = "pixel",
    labels = list(),
    filename = "synthetic_socint.csv",
    object.type = "TrackingData"
  )
}

test_that("approach is attributed to the animal that closed the distance", {
  socint <- load_socint()
  n <- 60
  fps <- 30
  # Animal 1 walks toward a stationary animal 2 at 30 px/s.
  a_x <- seq(0, by = 1, length.out = n)
  tracking <- socint_tracking(a_x, rep(0, n), rep(400, n), rep(0, n), fps = fps)

  config <- socint$config
  config$fps <- fps
  config$movement_cutoff <- 5
  config$approach_start_dist <- 1000
  config$retreat_start_dist <- 1000

  frames <- socint$make_event_table(tracking, config, "synthetic")

  expect_true(all(frames$a1_approaches_a2))
  expect_false(any(frames$a2_approaches_a1))
  expect_false(any(frames$a1_retreats_from_a2))
  expect_equal(unique(round(frames$a1_closing_speed, 6)), 30)
  expect_equal(unique(round(frames$a2_closing_speed, 6)), 0)
})

test_that("the animal that flees is credited with the separation", {
  socint <- load_socint()
  n <- 60
  fps <- 30
  # Animal 1 chases at 30 px/s; animal 2 outruns it at 60 px/s.
  a_x <- seq(0, by = 1, length.out = n)
  b_x <- seq(200, by = 2, length.out = n)
  tracking <- socint_tracking(a_x, rep(0, n), b_x, rep(0, n), fps = fps)

  config <- socint$config
  config$fps <- fps
  config$movement_cutoff <- 5
  config$approach_start_dist <- 1000
  config$retreat_start_dist <- 1000

  frames <- socint$make_event_table(tracking, config, "synthetic")

  # Animal 1 is approaching even though the pair separates.
  expect_true(all(frames$a1_approaches_a2))
  expect_true(all(frames$a2_retreats_from_a1))
  expect_true(all(frames$separation_led_by_a2))
  expect_false(any(frames$separation_led_by_a1))
  # The legacy speed-only rule fired both avoidance flags here because both
  # animals exceed the movement cutoff. The directed measures do not.
  expect_false(any(frames$a1_retreats_from_a2))
})

test_that("animals travelling in parallel produce no approach or retreat", {
  socint <- load_socint()
  n <- 60
  fps <- 30
  # Both animals move at 60 px/s on parallel tracks 100 px apart. Every
  # scalar speed cutoff is exceeded, but no relative motion occurs.
  step <- seq(0, by = 2, length.out = n)
  tracking <- socint_tracking(step, rep(0, n), step, rep(100, n), fps = fps)

  config <- socint$config
  config$fps <- fps
  config$movement_cutoff <- 5
  config$approach_start_dist <- 1000
  config$retreat_start_dist <- 1000
  config$avoidance_start_dist <- 1000

  frames <- socint$make_event_table(tracking, config, "synthetic")

  expect_equal(unique(round(frames$a1_closing_speed, 6)), 0)
  expect_equal(unique(round(frames$a2_closing_speed, 6)), 0)
  expect_false(any(frames$a1_approaches_a2))
  expect_false(any(frames$a2_approaches_a1))
  expect_false(any(frames$a1_retreats_from_a2))
  expect_false(any(frames$a2_retreats_from_a1))
  expect_false(any(frames$approach_event))
  expect_false(any(frames$retreat_event))
  expect_false(any(frames$a1_avoidance_from_a2_experimental))
  expect_false(any(frames$a2_avoidance_from_a1_experimental))
})

test_that("mutual approach is distinguished from one-sided approach", {
  socint <- load_socint()
  n <- 60
  fps <- 30
  a_x <- seq(0, by = 1, length.out = n)
  b_x <- seq(400, by = -1, length.out = n)
  tracking <- socint_tracking(a_x, rep(0, n), b_x, rep(0, n), fps = fps)

  config <- socint$config
  config$fps <- fps
  config$movement_cutoff <- 5
  config$approach_start_dist <- 1000
  config$retreat_start_dist <- 1000

  frames <- socint$make_event_table(tracking, config, "synthetic")

  expect_true(all(frames$both_approach))
  expect_true(all(frames$a1_approaches_a2))
  expect_true(all(frames$a2_approaches_a1))
  expect_equal(unique(round(frames$pair_closing_speed, 6)), 60)
})

test_that("the closing-speed decomposition sums to the pair closing speed", {
  socint <- load_socint()
  n <- 60
  fps <- 30
  a_x <- seq(0, by = 1.3, length.out = n)
  b_x <- seq(300, by = -0.4, length.out = n)
  tracking <- socint_tracking(a_x, rep(0, n), b_x, rep(0, n), fps = fps)

  config <- socint$config
  config$fps <- fps
  frames <- socint$make_event_table(tracking, config, "synthetic")

  expect_equal(
    frames$a1_closing_speed + frames$a2_closing_speed,
    frames$pair_closing_speed
  )
})

test_that("directed SocInt event columns are emitted with experimental labels", {
  socint <- load_socint()
  n <- 40
  tracking <- socint_tracking(
    seq(0, by = 1, length.out = n), rep(0, n),
    rep(300, n), rep(0, n), fps = 30
  )
  config <- socint$config
  config$fps <- 30
  frames <- socint$make_event_table(tracking, config, "synthetic")

  expect_true(all(c(
    "a1_approaches_a2", "a2_approaches_a1", "both_approach",
    "a1_retreats_from_a2", "a2_retreats_from_a1",
    "separation_led_by_a1", "separation_led_by_a2",
    "a1_closing_speed", "a2_closing_speed", "pair_closing_speed"
  ) %in% names(frames)))

  # The unvalidated heuristics must not be presented as established behaviors.
  expect_true(all(c(
    "a1_following_a2_experimental", "a2_following_a1_experimental",
    "a1_avoidance_from_a2_experimental", "a2_avoidance_from_a1_experimental"
  ) %in% names(frames)))
  expect_false("a1_following_a2" %in% names(frames))
  expect_false("a1_avoidance_from_a2" %in% names(frames))
})
