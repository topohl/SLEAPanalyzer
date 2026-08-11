#' @title DLCA_SocP v1.1.0.R
#' @description Batch analysis for the social-preference workflow.
#' @version 1.1.0 (Phase 1 correctness fixes)

required_packages <- c("sp", "imputeTS", "ggplot2", "cowplot", "stringr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required SocP package(s) before running: ", paste(missing_packages, collapse = ", "))
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
  batches = c("B3"),
  phases = c("HAB", "S1", "S2"),
  behavior_root = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior",
  animal_id_code_file = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Planning/animalIDCode.txt"
)

read_metadata_or_empty <- function(path, columns, sep = "") {
  if (!file.exists(path)) {
    warning("Metadata file not found: ", path)
    return(as.data.frame(setNames(replicate(length(columns), character(), simplify = FALSE), columns)))
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

animalIDCode <- read_metadata_or_empty(config$animal_id_code_file, c("Code", "ID"))

for (batch in config$batches) {
  for (socpPhase in config$phases) {
    inputDir <- file.path(config$behavior_root, batch, "SocP", "SLEAP", "formatted", socpPhase)
    outputDir <- file.path(config$behavior_root, batch, "SocP", "SLEAP", "output", socpPhase)
    novelLocPath <- file.path(config$behavior_root, batch, "SocP", paste0("novelLoc", socpPhase, ".txt"))
    plotDir <- file.path(outputDir, "plots")
    dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

    novelLoc <- read_metadata_or_empty(novelLocPath, c("Code", "NovelLoc"), sep = "\t")
    fileList <- list.files(path = inputDir, pattern = "\\.csv$", full.names = TRUE)
    dfList <- list()

    for (inputFile in fileList) {
      inputFileName <- tools::file_path_sans_ext(basename(inputFile))
      tracking <- ReadDLCDataFromCSV(file = inputFile, fps = config$fps)
      tracking <- CalibrateTrackingData(
        tracking, method = "area", in.metric = 44 * 24,
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
      metrics <- compute_socp_metrics(tracking, novel_location, config$fps)

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
          speedRaw = tracking$Report$bodycentre.raw.speed
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
        plot = plots$bodycentre, width = 7, height = 4
      )
      message("Processed file ", basename(inputFile))
    }

    if (length(dfList) == 0) {
      warning("No SocP CSV files found for ", batch, " ", socpPhase, " in ", inputDir)
      next
    }
    dfCombined <- do.call(rbind, dfList)
    utils::write.csv(dfCombined, file.path(outputDir, "combined_output.csv"), row.names = FALSE)
    message("Processing complete for ", batch, " ", socpPhase)
  }
}
