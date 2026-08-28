# Dependency-light geometry primitives shared by behavioral assays.

euclidean_distance <- function(x1, y1, x2, y2) {
  validate_numeric_vector(x1, "x1")
  validate_numeric_vector(y1, "y1")
  validate_numeric_vector(x2, "x2")
  validate_numeric_vector(y2, "y2")
  validate_equal_lengths(x1, y1, x2, y2, names = c("x1", "y1", "x2", "y2"))
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

vector_angle_degrees <- function(ax, ay, bx, by) {
  validate_numeric_vector(ax, "ax")
  validate_numeric_vector(ay, "ay")
  validate_numeric_vector(bx, "bx")
  validate_numeric_vector(by, "by")
  validate_equal_lengths(ax, ay, bx, by, names = c("ax", "ay", "bx", "by"))
  denominator <- sqrt(ax^2 + ay^2) * sqrt(bx^2 + by^2)
  cosine <- (ax * bx + ay * by) / denominator
  cosine[!is.finite(cosine)] <- NA_real_
  cosine <- pmax(-1, pmin(1, cosine))
  acos(cosine) * 180 / pi
}

tracking_point_distance <- function(tracking, first, second) {
  validate_tracking_data(tracking, required_landmarks = c(first, second))
  first_xy <- get_point_coordinates(tracking, first)
  second_xy <- get_point_coordinates(tracking, second)
  euclidean_distance(first_xy$x, first_xy$y, second_xy$x, second_xy$y)
}

tracking_vector_angle <- function(tracking, first_from, first_to,
                                  second_from, second_to) {
  required <- c(first_from, first_to, second_from, second_to)
  validate_tracking_data(tracking, required_landmarks = required)
  a_from <- get_point_coordinates(tracking, first_from)
  a_to <- get_point_coordinates(tracking, first_to)
  b_from <- get_point_coordinates(tracking, second_from)
  b_to <- get_point_coordinates(tracking, second_to)
  vector_angle_degrees(
    a_to$x - a_from$x,
    a_to$y - a_from$y,
    b_to$x - b_from$x,
    b_to$y - b_from$y
  )
}

points_in_axis_aligned_box <- function(x, y, center_x, center_y, width, height) {
  validate_numeric_vector(x, "x")
  validate_numeric_vector(y, "y")
  validate_numeric_vector(center_x, "center_x")
  validate_numeric_vector(center_y, "center_y")
  validate_equal_lengths(
    x, y, center_x, center_y,
    names = c("x", "y", "center_x", "center_y")
  )
  validate_scalar_number(width, "width", positive = TRUE)
  validate_scalar_number(height, "height", positive = TRUE)
  abs(x - center_x) <= width / 2 & abs(y - center_y) <= height / 2
}

validate_polygon <- function(polygon) {
  if (!is.data.frame(polygon) && !is.matrix(polygon)) {
    stop("polygon must be a data frame or matrix")
  }
  if (!all(c("x", "y") %in% colnames(polygon))) {
    stop("polygon must contain x and y columns")
  }
  if (nrow(polygon) < 3 || !is.numeric(polygon[, "x"]) ||
      !is.numeric(polygon[, "y"]) ||
      any(!is.finite(as.matrix(polygon[, c("x", "y"), drop = FALSE])))) {
    stop("polygon must contain at least three finite numeric vertices")
  }
  vertices <- data.frame(x = polygon[, "x"], y = polygon[, "y"])
  if (nrow(vertices) > 3 &&
      vertices$x[1] == vertices$x[nrow(vertices)] &&
      vertices$y[1] == vertices$y[nrow(vertices)]) {
    vertices <- vertices[-nrow(vertices), , drop = FALSE]
  }
  if (nrow(unique(vertices)) < 3) stop("polygon must contain at least three unique vertices")
  vertices
}

polygon_area <- function(polygon) {
  vertices <- validate_polygon(polygon)
  next_vertex <- c(2:nrow(vertices), 1L)
  abs(sum(vertices$x * vertices$y[next_vertex] -
            vertices$y * vertices$x[next_vertex])) / 2
}

points_in_polygon <- function(x, y, polygon, include_boundary = TRUE,
                              tolerance = sqrt(.Machine$double.eps)) {
  validate_numeric_vector(x, "x")
  validate_numeric_vector(y, "y")
  validate_equal_lengths(x, y, names = c("x", "y"))
  vertices <- validate_polygon(polygon)
  validate_scalar_number(tolerance, "tolerance", positive = TRUE, allow_zero = TRUE)

  inside <- rep(FALSE, length(x))
  boundary <- rep(FALSE, length(x))
  j <- nrow(vertices)
  for (i in seq_len(nrow(vertices))) {
    xi <- vertices$x[i]
    yi <- vertices$y[i]
    xj <- vertices$x[j]
    yj <- vertices$y[j]
    cross <- (x - xi) * (yj - yi) - (y - yi) * (xj - xi)
    on_segment <- abs(cross) <= tolerance &
      x >= min(xi, xj) - tolerance & x <= max(xi, xj) + tolerance &
      y >= min(yi, yj) - tolerance & y <= max(yi, yj) + tolerance
    boundary <- boundary | on_segment

    crosses <- ((yi > y) != (yj > y)) &
      (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
    crosses[is.na(crosses)] <- FALSE
    inside <- xor(inside, crosses)
    j <- i
  }
  if (include_boundary) inside | boundary else inside & !boundary
}

rectangular_arena_geometry <- function(corners) {
  vertices <- validate_polygon(corners)
  if (nrow(vertices) != 4) stop("rectangular arena geometry requires exactly four corners")
  next_vertex <- c(2:4, 1L)
  side_lengths <- euclidean_distance(
    vertices$x, vertices$y, vertices$x[next_vertex], vertices$y[next_vertex]
  )
  area <- polygon_area(vertices)
  list(
    valid = is.finite(area) && area > 0 && all(is.finite(side_lengths)) && all(side_lengths > 0),
    area = area,
    side_lengths = side_lengths,
    opposite_side_ratios = c(side_lengths[1] / side_lengths[3],
                             side_lengths[2] / side_lengths[4])
  )
}

#' Whether a polygon has any pair of non-adjacent edges that cross.
#'
#' Listing zone corners diagonally rather than around the perimeter produces a
#' self-intersecting shape. Point-in-polygon tests on such a shape are not
#' wrong in an obvious way; they quietly report most interior points as
#' outside, so occupancy is under-counted with no diagnostic.
#'
#' @param polygon a data frame or matrix with x and y columns
#' @return TRUE when any two non-adjacent edges intersect
is_self_intersecting <- function(polygon) {
  vertices <- validate_polygon(polygon)
  n <- nrow(vertices)
  if (n < 4) return(FALSE)

  orientation <- function(ax, ay, bx, by, cx, cy) {
    value <- (by - ay) * (cx - bx) - (bx - ax) * (cy - by)
    ifelse(abs(value) < 1e-12, 0, sign(value))
  }
  segments_cross <- function(p1, p2, p3, p4) {
    o1 <- orientation(p1[1], p1[2], p2[1], p2[2], p3[1], p3[2])
    o2 <- orientation(p1[1], p1[2], p2[1], p2[2], p4[1], p4[2])
    o3 <- orientation(p3[1], p3[2], p4[1], p4[2], p1[1], p1[2])
    o4 <- orientation(p3[1], p3[2], p4[1], p4[2], p2[1], p2[2])
    o1 != o2 && o3 != o4
  }

  at <- function(i) c(vertices$x[i], vertices$y[i])
  nxt <- function(i) if (i == n) 1L else i + 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (j <= i) next
      # Skip adjacent edges, which legitimately share a vertex.
      if (j == i || nxt(i) == j || nxt(j) == i) next
      if (segments_cross(at(i), at(nxt(i)), at(j), at(nxt(j)))) return(TRUE)
    }
  }
  FALSE
}
