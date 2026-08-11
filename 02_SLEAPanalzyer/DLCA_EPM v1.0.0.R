# Elevated-plus-maze batch workflow (Phase 1 correctness fixes).

required_packages <- c("sp", "imputeTS", "ggplot2", "cowplot")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required EPM package(s) before running: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    if (file.exists(script_path)) return(dirname(normalizePath(script_path)))
  }
  getwd()
}

script_dir <- get_script_dir()
source(file.path(script_dir, "DLCAnalyzer_Functions_final.R"))

config <- list(
  fps = 30,
  batches = c("B1", "B2", "B3", "B4", "B5", "B6"),
  behavior_root = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior",
  calibration_distance = 60,
  calibration_points = c("tl", "br"),
  zone_file = file.path(script_dir, "EPM_zoneinfo.csv")
)

if (!file.exists(config$zone_file)) stop("EPM zone file not found: ", config$zone_file)
zoneInfo <- utils::read.table(
  config$zone_file, sep = ";", header = TRUE,
  stringsAsFactors = FALSE, check.names = FALSE
)
zone_points <- unique(unlist(zoneInfo, use.names = FALSE))
zone_points <- zone_points[!is.na(zone_points) & nzchar(trimws(zone_points))]

for (batch in config$batches) {
  inputFolder <- file.path(config$behavior_root, batch, "EPM", "SLEAP", "formatted")
  outputFolder <- file.path(config$behavior_root, batch, "EPM", "SLEAP", "output")
  overviewPlotFolder <- file.path(config$behavior_root, batch, "EPM", "SLEAP", "OverviewPlot")
  dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
  dir.create(overviewPlotFolder, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(inputFolder, pattern = "\\.csv$", full.names = TRUE)
  pipeline <- function(path) {
    tracking <- ReadDLCDataFromCSV(file = path, fps = config$fps)
    required_points <- unique(c(
      config$calibration_points, zone_points,
      "headcentre", "bodycentre", "neck"
    ))
    missing_points <- setdiff(required_points, names(tracking$data))
    if (length(missing_points) > 0) {
      stop(basename(path), " is missing EPM point(s): ", paste(missing_points, collapse = ", "))
    }
    tracking <- CalibrateTrackingData(
      tracking, method = "distance", in.metric = config$calibration_distance,
      points = config$calibration_points
    )
    tracking <- AddZones(tracking, zoneInfo)
    EPMAnalysis(
      tracking, movement_cutoff = 5, integration_period = 5,
      points = "bodycentre", nosedips = TRUE
    )
  }

  if (length(files) == 0) {
    warning("No EPM CSV files found for batch ", batch, " in ", inputFolder)
    next
  }

  trackingAll <- lapply(files, pipeline)
  names(trackingAll) <- basename(files)
  report <- MultiFileReport(trackingAll)
  utils::write.csv(report, file.path(outputFolder, "Report.csv"), row.names = FALSE)

  for (i in seq_along(files)) {
    plot <- OverviewPlot(trackingAll[[i]], "bodycentre")
    outputFile <- file.path(
      overviewPlotFolder,
      paste0(tools::file_path_sans_ext(basename(files[i])), ".tiff")
    )
    ggplot2::ggsave(
      filename = outputFile, plot = plot, device = "tiff",
      width = 8, height = 10
    )
  }
  message("EPM processing complete for batch ", batch)
}
