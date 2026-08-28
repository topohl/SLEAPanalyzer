# Synthetic dyads whose directional attribution is known analytically.

test_that("point_velocity recovers a known constant velocity", {
  n <- 11
  # 2 units per frame in x at 10 fps is 20 units per second.
  v <- point_velocity(seq(0, by = 2, length.out = n), rep(0, n), fps = 10)
  expect_equal(v$vx, rep(20, n))
  expect_equal(v$vy, rep(0, n))
})

test_that("only A moves: the approach is attributed entirely to A", {
  n <- 11
  ax <- seq(0, by = 1, length.out = n)   # A walks toward B at 10 units/s
  motion <- dyadic_relative_motion(ax, rep(0, n), rep(100, n), rep(0, n), fps = 10)

  expect_equal(motion$a_closing, rep(10, n))
  expect_equal(motion$b_closing, rep(0, n))
  expect_equal(motion$closing_speed, rep(10, n))
  # The decomposition matches the observed change in distance.
  expect_equal(motion$closing_speed, motion$observed_closing_speed)

  labels <- classify_relative_motion(motion, threshold = 1)
  expect_true(all(labels$a_approaches_b))
  expect_false(any(labels$b_approaches_a))
  expect_true(all(labels$approach_led_by_a))
  expect_false(any(labels$approach_led_by_b))
})

test_that("only B moves: the approach is attributed entirely to B", {
  n <- 11
  bx <- seq(100, by = -1, length.out = n)
  motion <- dyadic_relative_motion(rep(0, n), rep(0, n), bx, rep(0, n), fps = 10)

  expect_equal(motion$a_closing, rep(0, n))
  expect_equal(motion$b_closing, rep(10, n))

  labels <- classify_relative_motion(motion, threshold = 1)
  expect_false(any(labels$a_approaches_b))
  expect_true(all(labels$b_approaches_a))
  expect_true(all(labels$approach_led_by_b))
})

test_that("both animals approach each other", {
  n <- 11
  ax <- seq(0, by = 1, length.out = n)
  bx <- seq(100, by = -1, length.out = n)
  motion <- dyadic_relative_motion(ax, rep(0, n), bx, rep(0, n), fps = 10)

  expect_equal(motion$a_closing, rep(10, n))
  expect_equal(motion$b_closing, rep(10, n))
  expect_equal(motion$closing_speed, rep(20, n))

  labels <- classify_relative_motion(motion, threshold = 1)
  expect_true(all(labels$both_approach))
  expect_true(all(labels$pair_approaching))
  expect_false(any(labels$pair_separating))
})

test_that("A pursues while B flees: attribution separates chaser from fleer", {
  n <- 11
  # Both animals move right. A starts behind at 0 and moves at 10 units/s;
  # B starts ahead at 100 and outruns it at 20 units/s, so the pair separates.
  ax <- seq(0, by = 1, length.out = n)
  bx <- seq(100, by = 2, length.out = n)
  motion <- dyadic_relative_motion(ax, rep(0, n), bx, rep(0, n), fps = 10)

  # A is closing on B; B is opening the gap. The pair separates net.
  expect_equal(motion$a_closing, rep(10, n))
  expect_equal(motion$b_closing, rep(-20, n))
  expect_equal(motion$closing_speed, rep(-10, n))
  expect_equal(motion$closing_speed, motion$observed_closing_speed)

  labels <- classify_relative_motion(motion, threshold = 1)
  expect_true(all(labels$a_approaches_b))
  expect_true(all(labels$b_retreats_from_a))
  expect_false(any(labels$a_retreats_from_b))
  expect_false(any(labels$b_approaches_a))
  expect_true(all(labels$pair_separating))
  # The decisive property: the separation is attributed to B, the animal that
  # actually produced it. Both animals exceed any scalar speed cutoff, so the
  # legacy speed-only rule would have flagged both as avoiding.
  expect_true(all(labels$separation_led_by_b))
  expect_false(any(labels$separation_led_by_a))
})

test_that("separation driven by one animal is attributed to that animal", {
  n <- 11
  # A flees to the left at 20 units/s; B drifts slowly toward A at 5 units/s.
  ax <- seq(0, by = -2, length.out = n)
  bx <- seq(100, by = -0.5, length.out = n)
  motion <- dyadic_relative_motion(ax, rep(0, n), bx, rep(0, n), fps = 10)

  expect_equal(motion$a_closing, rep(-20, n))
  expect_equal(motion$b_closing, rep(5, n))
  expect_equal(motion$closing_speed, rep(-15, n))

  labels <- classify_relative_motion(motion, threshold = 1)
  expect_true(all(labels$pair_separating))
  expect_true(all(labels$separation_led_by_a))
  expect_false(any(labels$separation_led_by_b))
  expect_true(all(labels$a_retreats_from_b))
  # B is still approaching even though the pair separates.
  expect_true(all(labels$b_approaches_a))
})

test_that("both animals moving in parallel produce no relative motion", {
  n <- 11
  # A speed-only rule would call both animals "moving" and fire both
  # avoidance flags; the directional decomposition correctly reports nothing.
  ax <- seq(0, by = 2, length.out = n)
  bx <- seq(0, by = 2, length.out = n)
  motion <- dyadic_relative_motion(ax, rep(0, n), bx, rep(30, n), fps = 10)

  expect_equal(motion$a_closing, rep(0, n))
  expect_equal(motion$b_closing, rep(0, n))
  expect_equal(motion$closing_speed, rep(0, n))
  expect_equal(motion$a_speed, rep(20, n))
  expect_equal(motion$b_speed, rep(20, n))

  labels <- classify_relative_motion(motion, threshold = 1)
  expect_false(any(labels$a_approaches_b))
  expect_false(any(labels$b_approaches_a))
  expect_false(any(labels$a_retreats_from_b))
  expect_false(any(labels$b_retreats_from_a))
  expect_false(any(labels$pair_approaching))
  expect_false(any(labels$pair_separating))
})

test_that("motion perpendicular to the axis does not change distance", {
  n <- 11
  motion <- dyadic_relative_motion(
    rep(0, n), seq(0, by = 1, length.out = n),
    rep(50, n), rep(0, n),
    fps = 10
  )
  # A moves at right angles to the line joining them at the first frame.
  expect_equal(motion$a_closing[1], 0)
  expect_equal(motion$b_closing[1], 0)
})

test_that("the decomposition matches the observed closing speed on straight paths", {
  n <- 21
  ax <- seq(0, by = 1.5, length.out = n)
  bx <- seq(200, by = -0.7, length.out = n)
  motion <- dyadic_relative_motion(ax, rep(0, n), bx, rep(0, n), fps = 25)
  expect_equal(motion$closing_speed, motion$observed_closing_speed, tolerance = 1e-8)
})

test_that("coincident animals yield NA rather than a fabricated direction", {
  n <- 5
  motion <- dyadic_relative_motion(rep(0, n), rep(0, n), rep(0, n), rep(0, n), fps = 10)
  expect_true(all(is.na(motion$a_closing)))
  expect_true(all(is.na(motion$b_closing)))

  labels <- classify_relative_motion(motion, threshold = 1)
  expect_true(all(is.na(labels$a_approaches_b)))
})

test_that("the threshold suppresses sub-threshold drift", {
  n <- 11
  ax <- seq(0, by = 0.01, length.out = n)
  motion <- dyadic_relative_motion(ax, rep(0, n), rep(100, n), rep(0, n), fps = 10)
  expect_equal(motion$a_closing, rep(0.1, n))

  labels <- classify_relative_motion(motion, threshold = 1)
  expect_false(any(labels$a_approaches_b))
  expect_true(all(classify_relative_motion(motion, threshold = 0.05)$a_approaches_b))
})

test_that("classify_relative_motion validates its inputs", {
  motion <- dyadic_relative_motion(0:4, rep(0, 5), rep(10, 5), rep(0, 5), fps = 10)
  expect_error(classify_relative_motion(motion, threshold = -1), "threshold")
  expect_error(classify_relative_motion(data.frame(a = 1), threshold = 1), "missing column")
})
