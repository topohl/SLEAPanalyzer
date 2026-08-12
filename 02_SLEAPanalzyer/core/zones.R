# Generic zone membership built on the shared polygon primitives.

point_in_zone <- function(x, y, zone, include_boundary = TRUE) {
  points_in_polygon(x, y, zone, include_boundary = include_boundary)
}

tracking_point_in_zone <- function(tracking, point, zone, include_boundary = TRUE) {
  coordinates <- get_point_coordinates(tracking, point)
  point_in_zone(coordinates$x, coordinates$y, zone, include_boundary = include_boundary)
}

tracking_point_in_named_zone <- function(tracking, point, zone_name,
                                         include_boundary = TRUE) {
  if (is.null(tracking$zones) || !zone_name %in% names(tracking$zones)) {
    stop("TrackingData does not contain zone: ", zone_name)
  }
  tracking_point_in_zone(
    tracking, point, tracking$zones[[zone_name]], include_boundary = include_boundary
  )
}
