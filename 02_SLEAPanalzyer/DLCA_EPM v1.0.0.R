#' @title DLCA_EPM v1.0.0.R
#' @description Batch analysis for the elevated-plus-maze workflow.
#'
#' Configuration is read from a YAML file, not from this script:
#'
#'   SLEAP_ANALYZER_CONFIG=my_epm.yaml Rscript "DLCA_EPM v1.0.0.R"
#'
#' See config/epm.example.yaml for a documented template.

required_packages <- c("sp", "ggplot2", "cowplot", "yaml")
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
source(file.path(script_dir, "core", "assay_config.R"))

config <- load_assay_config(resolve_config_path(assay = "EPM", script_dir = script_dir), "EPM")
repo_dir <- dirname(script_dir)

if (!file.exists(config$zone_file)) stop("EPM zone file not found: ", config$zone_file)
zoneInfo <- utils::read.table(
  config$zone_file, sep = ";", header = TRUE,
  stringsAsFactors = FALSE, check.names = FALSE
)
zone_points <- unique(unlist(zoneInfo, use.names = FALSE))
zone_points <- zone_points[!is.na(zone_points) & nzchar(trimws(zone_points))]

# Resolve the calibration once, so a misconfigured run fails before any file
# is read rather than part way through a batch.
#
# `distance` is the documented EPM default because the maze is a plus, not a
# rectangle: see the note in core/assay_config.R. `area` stays available for
# configurations that deliberately supply a polygon enclosing a known area,
# in which case `calibration_points` must trace that polygon in order around
# its perimeter -- listing them diagonally yields a self-intersecting shape
# whose area is meaningless.
calibration <- list(
  method = if (is.null(config$calibration_method)) "distance" else config$calibration_method
)
if (calibration$method == "distance") {
  calibration$points <- if (is.null(config$calibration_points)) {
    stop("calibration_method = 'distance' requires calibration_points (two landmarks)")
  } else config$calibration_points
  if (length(calibration$points) != 2) {
    stop("calibration_method = 'distance' requires exactly two calibration_points, got ",
         length(calibration$points))
  }
  if (is.null(config$calibration_distance_cm)) {
    stop("calibration_method = 'distance' requires calibration_distance_cm")
  }
  calibration$in.metric <- config$calibration_distance_cm
} else {
  calibration$points <- if (is.null(config$calibration_points)) {
    config$arena_corner_names
  } else config$calibration_points
  if (length(calibration$points) < 3) {
    stop("calibration_method = 'area' requires at least three calibration_points")
  }
  calibration$in.metric <- config$arena_width_cm * config$arena_height_cm
}

batches <- if (is.null(config$batches)) "" else config$batches

for (batch in batches) {
  inputFolder <- if (nzchar(batch)) file.path(config$input_dir, batch) else config$input_dir
  outputFolder <- if (nzchar(batch)) file.path(config$output_dir, batch) else config$output_dir
  overviewPlotFolder <- file.path(outputFolder, "OverviewPlot")
  dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
  dir.create(overviewPlotFolder, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(inputFolder, pattern = "\\.csv$", full.names = TRUE)
  qcList <- list()

  pipeline <- function(path) {
    tracking <- read_tracking_csv(file = path, fps = config$fps)
    required_points <- unique(c(
      calibration$points, zone_points,
      "headcentre", "bodycentre", "neck"
    ))
    missing_points <- setdiff(required_points, names(tracking$data))
    if (length(missing_points) > 0) {
      stop(basename(path), " is missing EPM point(s): ", paste(missing_points, collapse = ", "))
    }

    # Bounded interpolation before analysis: EPM previously ran on raw
    # coordinates, so untracked frames were scored as being outside every
    # zone and, through the inverted-zone path, sometimes inside one.
    tracking <- interpolate_tracking(
      tracking,
      landmarks = c("headcentre", "bodycentre", "neck"),
      max_gap_s = config$max_interpolation_gap_s,
      likelihood_cutoff = config$likelihood_cutoff
    )

    tracking <- CalibrateTrackingData(
      tracking, method = calibration$method,
      in.metric = calibration$in.metric,
      points = calibration$points
    )
    tracking <- AddZones(tracking, zoneInfo)
    EPMAnalysis(
      tracking,
      movement_cutoff = config$movement_cutoff_cm_s,
      integration_period = config$integration_period_frames,
      points = "bodycentre",
      nosedips = isTRUE(config$nose_dips),
      nosedip_integration_period = if (is.null(config$nosedip_integration_period_frames)) {
        config$integration_period_frames
      } else {
        config$nosedip_integration_period_frames
      }
    )
  }

  if (length(files) == 0) {
    warning("No EPM CSV files found for batch ", batch, " in ", inputFolder)
    next
  }

  trackingAll <- lapply(files, pipeline)
  names(trackingAll) <- basename(files)

  for (tracking in trackingAll) {
    qcReport <- tracking_qc_report(
      tracking,
      landmarks = c("headcentre", "bodycentre", "neck"),
      required_landmarks = "bodycentre",
      likelihood_cutoff = config$likelihood_cutoff,
      max_speed = config$max_plausible_speed_cm_s
    )
    qcDecision <- qc_flags(
      qcReport,
      min_valid_fraction = config$qc_min_valid_fraction,
      max_longest_gap_s = config$qc_max_longest_gap_s,
      max_interpolated_fraction = config$qc_max_interpolated_fraction
    )
    if (!qcDecision$pass) {
      warning(
        tracking$filename, " failed tracking QC: ",
        paste(qcDecision$reasons, collapse = "; "),
        ". The result is still written; exclusion is an explicit decision."
      )
    }
    qcList[[length(qcList) + 1L]] <- cbind(
      qc_summary_row(qcReport),
      data.frame(
        qcPass = qcDecision$pass,
        qcReasons = paste(qcDecision$reasons, collapse = "; "),
        stringsAsFactors = FALSE
      )
    )
  }

  report <- MultiFileReport(trackingAll)
  utils::write.csv(report, file.path(outputFolder, "Report.csv"), row.names = FALSE)
  utils::write.csv(
    do.call(rbind, qcList), file.path(outputFolder, "tracking_qc.csv"), row.names = FALSE
  )

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

  if (isTRUE(config$write_manifest)) {
    manifest <- run_manifest(
      config = config,
      inputs = files,
      packages = required_packages,
      repo_dir = repo_dir,
      extra = list(assay = "EPM", batch = batch, files = length(files))
    )
    assert_reproducible_run(manifest)
    write_run_manifest(manifest, file.path(outputFolder, "run_manifest.yaml"))
  }
  message("EPM processing complete for batch ", batch)
}
