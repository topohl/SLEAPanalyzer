# Canonical dyadic geometry and directional relative motion.
#
# The problem this module solves: scalar speed cannot attribute a change in
# inter-animal distance to an animal. If two animals separate and both are
# moving, a speed-only rule marks both as avoiding, which is not a measurement.
#
# The decomposition used here projects each animal's velocity onto the
# inter-animal axis, which splits the observed closing speed into the part
# each animal contributed.

#' Per-frame velocity of a tracked point.
#'
#' Central differences are used for interior frames and one-sided differences
#' at the ends, so the result is in coordinate units per second and is not
#' shifted by half a frame relative to the position series.
#'
#' @param x,y coordinate vectors
#' @param fps frames per second
#' @return a list of vx and vy in coordinate units per second
point_velocity <- function(x, y, fps) {
  validate_numeric_vector(x, "x")
  validate_numeric_vector(y, "y")
  validate_equal_lengths(x, y, names = c("x", "y"))
  validate_fps(fps)
  n <- length(x)
  if (n < 2) {
    return(list(vx = rep(NA_real_, n), vy = rep(NA_real_, n)))
  }
  central <- function(v) {
    out <- rep(NA_real_, n)
    out[1] <- (v[2] - v[1]) * fps
    out[n] <- (v[n] - v[n - 1]) * fps
    if (n > 2) {
      interior <- 2:(n - 1)
      out[interior] <- (v[interior + 1] - v[interior - 1]) * fps / 2
    }
    out
  }
  list(vx = central(x), vy = central(y))
}

#' Directional relative motion between two animals.
#'
#' Let u be the unit vector from animal A to animal B and d the distance
#' between them. Then
#'
#'   d(d)/dt = (v_B - v_A) . u
#'
#' so the closing speed (positive when the animals approach) decomposes as
#'
#'   closing = v_A . u  +  (-v_B . u)
#'           = a_closing + b_closing
#'
#' `a_closing` is the part of the approach produced by animal A moving toward
#' animal B, and `b_closing` the part produced by animal B moving toward
#' animal A. Either term can be negative, which means that animal is moving
#' away. This is what makes "A approached B" separable from "B approached A".
#'
#' In discrete time the decomposed and observed closing speeds differ by a
#' curvature term because u is evaluated at the current frame; they agree for
#' straight-line motion and are reported separately so the discrepancy is
#' visible rather than hidden.
#'
#' @param ax,ay animal A coordinates
#' @param bx,by animal B coordinates
#' @param fps frames per second
#' @return a data frame with one row per frame
dyadic_relative_motion <- function(ax, ay, bx, by, fps) {
  validate_equal_lengths(ax, ay, bx, by, names = c("ax", "ay", "bx", "by"))
  validate_fps(fps)

  distance <- euclidean_distance(ax, ay, bx, by)
  # Unit vector from A to B. Undefined when the points coincide.
  ux <- (bx - ax) / distance
  uy <- (by - ay) / distance
  ux[!is.finite(ux)] <- NA_real_
  uy[!is.finite(uy)] <- NA_real_

  va <- point_velocity(ax, ay, fps)
  vb <- point_velocity(bx, by, fps)

  a_closing <- va$vx * ux + va$vy * uy
  b_closing <- -(vb$vx * ux + vb$vy * uy)

  observed_closing <- rep(NA_real_, length(distance))
  if (length(distance) >= 2) {
    n <- length(distance)
    observed_closing[1] <- -(distance[2] - distance[1]) * fps
    observed_closing[n] <- -(distance[n] - distance[n - 1]) * fps
    if (n > 2) {
      interior <- 2:(n - 1)
      observed_closing[interior] <- -(distance[interior + 1] - distance[interior - 1]) * fps / 2
    }
  }

  data.frame(
    distance = distance,
    a_closing = a_closing,
    b_closing = b_closing,
    closing_speed = a_closing + b_closing,
    observed_closing_speed = observed_closing,
    a_speed = sqrt(va$vx^2 + va$vy^2),
    b_speed = sqrt(vb$vx^2 + vb$vy^2)
  )
}

#' Attribute approach and retreat to individual animals.
#'
#' Each animal is classified independently from its own contribution to the
#' closing speed, so "both approach", "A approaches while B retreats" and
#' "separation driven mainly by B" are all distinguishable.
#'
#' `threshold` is a speed in coordinate units per second and must match the
#' coordinate unit of the tracking data.
#'
#' @param motion the result of dyadic_relative_motion()
#' @param threshold minimum directional speed counted as approach or retreat
#' @return a data frame of logical per-frame classifications
classify_relative_motion <- function(motion, threshold) {
  validate_scalar_number(threshold, "threshold", positive = TRUE, allow_zero = TRUE)
  required <- c("a_closing", "b_closing", "closing_speed")
  missing <- setdiff(required, names(motion))
  if (length(missing) > 0) {
    stop("motion is missing column(s): ", paste(missing, collapse = ", "))
  }

  a_approaches <- motion$a_closing > threshold
  b_approaches <- motion$b_closing > threshold
  a_retreats <- motion$a_closing < -threshold
  b_retreats <- motion$b_closing < -threshold
  approaching <- motion$closing_speed > threshold
  separating <- motion$closing_speed < -threshold

  # Which animal contributed more to the change in distance. Only meaningful
  # while the pair is actually approaching or separating.
  a_dominant <- abs(motion$a_closing) > abs(motion$b_closing)

  data.frame(
    a_approaches_b = a_approaches,
    b_approaches_a = b_approaches,
    a_retreats_from_b = a_retreats,
    b_retreats_from_a = b_retreats,
    both_approach = a_approaches & b_approaches,
    both_retreat = a_retreats & b_retreats,
    pair_approaching = approaching,
    pair_separating = separating,
    approach_led_by_a = approaching & a_dominant,
    approach_led_by_b = approaching & !a_dominant,
    separation_led_by_a = separating & a_dominant,
    separation_led_by_b = separating & !a_dominant
  )
}
