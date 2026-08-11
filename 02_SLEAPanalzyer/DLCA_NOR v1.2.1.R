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

read_metadata_or_empty <- function(path, columns, sep = "") {
  if (!file.exists(path)) {
    warning("Metadata file not found: ", path)
    out <- as.data.frame(setNames(replicate(length(columns), character(), simplify = FALSE), columns))
    return(out)
  }
  out <- utils::read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  missing <- setdiff(columns, names(out))
  if (length(missing) > 0) stop("Metadata file is missing column(s): ", paste(missing, collapse = ", "))
  out
}

metadata_value <- function(data, code, column) {
  values <- data[data$Code == code, column]
  if (length(values) == 0) return(NA_character_)
  if (length(values) > 1) warning("Multiple metadata rows found for code ", code, "; using the first")
  as.character(values[1])
}

fill_edges_and_gaps <- function(x) {
  x <- zoo::na.locf(x, na.rm = FALSE)
  zoo::na.locf(x, fromLast = TRUE, na.rm = FALSE)
}

animalIDCode <- read_metadata_or_empty(config$animal_id_code_file, c("Code", "ID"))

for (batch in config$batches) {
  inputDir <- file.path(config$behavior_root, batch, "NOR", "SLEAP", "formatted")
  outputDir <- file.path(config$behavior_root, batch, "NOR", "SLEAP", "output")
  novelLocPath <- file.path(config$behavior_root, batch, "NOR", config$novel_location_file)
  plotDir <- file.path(outputDir, "plots")
  dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

  novelLoc <- read_metadata_or_empty(novelLocPath, c("Code", "NovelLoc"), sep = "\t")
  fileList <- list.files(path = inputDir, pattern = "\\.csv$", full.names = TRUE)
  dfList <- list()

  for (inputFile in fileList) {
    inputFileName <- tools::file_path_sans_ext(basename(inputFile))
    tracking <- ReadDLCDataFromCSV(file = inputFile, fps = config$fps)

    for (point in c("nose", "bodycentre")) {
      if (!point %in% names(tracking$data)) stop(inputFileName, " is missing tracked point: ", point)
      tracking$data[[point]]$x <- fill_edges_and_gaps(tracking$data[[point]]$x)
      tracking$data[[point]]$y <- fill_edges_and_gaps(tracking$data[[point]]$y)
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
    novel_location <- metadata_value(novelLoc, code, "NovelLoc")
    metrics <- compute_nor_metrics(tracking, novel_location, config$fps)

    spine1_distance <- if (all(c("spine1", "bodycentre") %in% names(tracking$data))) {
      tracking_distance(tracking, "spine1", "bodycentre")
    } else numeric()
    spine2_distance <- if (all(c("spine2", "bodycentre") %in% names(tracking$data))) {
      tracking_distance(tracking, "bodycentre", "spine2")
    } else numeric()
    frequencyRear <- if (length(spine1_distance) == 0 || length(spine2_distance) == 0) {
      NA_integer_
    } else {
      event_entry_count(spine1_distance <= 1 & spine2_distance <= 1)
    }

    df <- cbind(
      data.frame(
        file = inputFileName,
        ID = metadata_value(animalIDCode, code, "ID"),
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
