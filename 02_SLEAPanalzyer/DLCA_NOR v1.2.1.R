#' @title DLCA_NOR v1.2.1.R
#' @description Batch analysis for the novel-object-recognition workflow.
#' @version 1.2.1 (Phase 1 correctness fixes)

required_packages <- c("sp", "imputeTS", "ggplot2", "cowplot", "zoo", "stringr", "openxlsx")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required NOR package(s) before running: ", paste(missing_packages, collapse = ", "))
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
source(file.path(script_dir, "Behavioral_Metrics_Phase1.R"))

config <- list(
  fps = 30,
  batches = c("B1", "B2", "B3", "B4", "B5", "B6"),
  behavior_root = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior",
  animal_id_code_file = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Planning/animalIDCode.txt",
  novel_location_file = "novelLoc.txt"
)

fill_edges_and_gaps <- function(x) {
  x <- zoo::na.locf(x, na.rm = FALSE)
  zoo::na.locf(x, fromLast = TRUE, na.rm = FALSE)
}

animalIDCode <- read_metadata_table(
  config$animal_id_code_file, c("Code", "ID"), if_missing = "empty"
)

for (batch in config$batches) {
  inputDir <- file.path(config$behavior_root, batch, "NOR", "SLEAP", "formatted")
  outputDir <- file.path(config$behavior_root, batch, "NOR", "SLEAP", "output")
  novelLocPath <- file.path(config$behavior_root, batch, "NOR", config$novel_location_file)
  plotDir <- file.path(outputDir, "plots")
  dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

  novelLoc <- read_metadata_table(
    novelLocPath, c("Code", "NovelLoc"), sep = "\t", if_missing = "empty"
  )
  fileList <- list.files(path = inputDir, pattern = "\\.csv$", full.names = TRUE)
  dfList <- list()

  for (inputFile in fileList) {
    inputFileName <- tools::file_path_sans_ext(basename(inputFile))
    tracking <- read_tracking_csv(file = inputFile, fps = config$fps)

    for (point in c("nose", "bodycentre")) {
      if (!has_landmarks(tracking, point)) stop(inputFileName, " is missing tracked point: ", point)
      coordinates <- get_point_coordinates(tracking, point)
      coordinates$x <- fill_edges_and_gaps(coordinates$x)
      coordinates$y <- fill_edges_and_gaps(coordinates$y)
      tracking <- set_point_coordinates(tracking, point, coordinates)
    }

    tracking <- CalibrateTrackingData(
      tracking, method = "area", in.metric = 49 * 49,
      points = c("tl", "tr", "br", "bl")
    )
    tracking <- AddOFTZones(
      tracking, scale_center = 0.5, scale_periphery = 0.8,
      scale_corners = 0.4, points = c("tl", "tr", "br", "bl")
    )
    tracking <- OFTAnalysis(
      tracking, points = "bodycentre", movement_cutoff = 5,
      integration_period = 5
    )

    code <- stringr::str_extract(inputFileName, "^[A-Za-z0-9]{4}")
    if (is.na(code)) warning("Could not extract a four-character animal code from ", inputFileName)
    novel_location <- metadata_lookup(novelLoc, code, "NovelLoc")
    metrics <- compute_nor_metrics(tracking, novel_location, config$fps)

    spine1_distance <- if (has_landmarks(tracking, c("spine1", "bodycentre"))) {
      tracking_point_distance(tracking, "spine1", "bodycentre")
    } else numeric()
    spine2_distance <- if (has_landmarks(tracking, c("spine2", "bodycentre"))) {
      tracking_point_distance(tracking, "bodycentre", "spine2")
    } else numeric()
    frequencyRear <- if (length(spine1_distance) == 0 || length(spine2_distance) == 0) {
      NA_integer_
    } else {
      event_frequency(spine1_distance <= 1 & spine2_distance <= 1)
    }

    df <- cbind(
      data.frame(
        file = inputFileName,
        ID = metadata_lookup(animalIDCode, code, "ID"),
        Code = code,
        stringsAsFactors = FALSE
      ),
      metrics$summary,
      data.frame(
        distance = tracking$Report$bodycentre.raw.distance,
        stationary = tracking$Report$bodycentre.time.stationary,
        speedMoving = tracking$Report$bodycentre.speed.moving,
        speedRaw = tracking$Report$bodycentre.raw.speed,
        frequencyRear = frequencyRear
      )
    )

    utils::write.csv(
      df, file.path(outputDir, paste0(inputFileName, "_output.csv")),
      row.names = FALSE
    )
    dfList[[length(dfList) + 1L]] <- df

    plots <- PlotDensityPaths(tracking, points = "bodycentre")
    ggplot2::ggsave(
      file.path(plotDir, paste0(inputFileName, "_DensityPath.png")),
      plot = plots$bodycentre, width = 7, height = 6
    )
    message("Processed file ", basename(inputFile))
  }

  if (length(dfList) == 0) {
    warning("No NOR CSV files found for batch ", batch, " in ", inputDir)
    next
  }
  dfCombined <- do.call(rbind, dfList)
  openxlsx::write.xlsx(dfCombined, file.path(outputDir, "combined_output.xlsx"), rowNames = FALSE)
  message("Processing complete for NOR batch ", batch)
}
