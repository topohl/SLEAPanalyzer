repo_root <- Sys.getenv("SLEAP_ANALYZER_REPO_ROOT")
if (!nzchar(repo_root)) stop("SLEAP_ANALYZER_REPO_ROOT is not set; use tests/testthat.R")
source(file.path(repo_root, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(repo_root, "02_SLEAPanalzyer", "Behavioral_Metrics_Phase1.R"))

make_tracking <- function(frames = 0:9, fps = 30, points = list()) {
  n <- length(frames)
  if (length(points) == 0) points <- list(bodycentre = cbind(x = rep(0, n), y = rep(0, n)))
  data <- lapply(points, function(xy) {
    stopifnot(nrow(xy) == n)
    data.frame(
      frame = frames,
      x = as.numeric(xy[, "x"]),
      y = as.numeric(xy[, "y"]),
      likelihood = rep(1, n)
    )
  })
  median_data <- do.call(rbind, lapply(names(data), function(point) {
    data.frame(
      PointName = point,
      x = median(data[[point]]$x, na.rm = TRUE),
      y = median(data[[point]]$y, na.rm = TRUE),
      row.names = point
    )
  }))
  list(
    data = data,
    frames = frames,
    fps = fps,
    seconds = (frames - frames[1]) / fps,
    median.data = median_data,
    point.info = data.frame(PointName = names(data), PointType = "NotDefined"),
    distance.units = "pixel",
    labels = list(),
    filename = "synthetic.csv",
    object.type = "TrackingData"
  )
}

square_points <- function(n = 10, scale = 1, offset = c(0, 0)) {
  at <- function(x, y) cbind(
    x = rep(offset[1] + scale * x, n),
    y = rep(offset[2] + scale * y, n)
  )
  list(
    tl = at(0, 100), tr = at(100, 100),
    br = at(100, 0), bl = at(0, 0)
  )
}
