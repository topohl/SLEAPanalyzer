#' @title DLCA_SocP v1.1.0.R
#' @description Batch analysis for the social-preference workflow.
#' @version 1.1.0 (Phase 1 correctness fixes)

required_packages <- c("sp", "ggplot2", "stringr")
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
  animal_id_code_file = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Planning/animalIDCode.txt",

  # Arena geometry of the three-chamber apparatus floor.
  arena_width_cm = 44,
  arena_height_cm = 24,
  arena_corner_names = c("tl", "tr", "br", "bl"),

  # Tracking quality.
  max_interpolation_gap_s = 0.2,
  likelihood_cutoff = NULL,

  # Contact definition. Requires assay-specific validation against manually
  # scored video; see docs/assay_definitions.md.
  contact_distance_cm = 6,
  proximity_range_cm = c(6, 10),
  min_contact_bout_s = 0,
  max_contact_gap_s = 0
)

# Metadata parsing uses the shared core helpers; SocP previously carried
# private copies that could drift from the NOR implementation.
animalIDCode <- read_metadata_table(
  config$animal_id_code_file, c("Code", "ID"), if_missing = "empty"
)

for (batch in config$batches) {
  for (socpPhase in config$phases) {
    inputDir <- file.path(config$behavior_root, batch, "SocP", "SLEAP", "formatted", socpPhase)
    outputDir <- file.path(config$behavior_root, batch, "SocP", "SLEAP", "output", socpPhase)
    novelLocPath <- file.path(config$behavior_root, batch, "SocP", paste0("novelLoc", socpPhase, ".txt"))
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

      # Bounded interpolation: SocP previously performed none at all, so
      # untracked frames reached the contact test as NA and were scored as
      # confident absence of social contact.
      tracking <- interpolate_tracking(
        tracking,
        landmarks = c("nose", "bodycentre"),
        max_gap_s = config$max_interpolation_gap_s,
        likelihood_cutoff = config$likelihood_cutoff
      )
      trackingQC <- interpolation_report(tracking, c("nose", "bodycentre"))

      tracking <- CalibrateTrackingData(
        tracking, method = "area",
        in.metric = config$arena_width_cm * config$arena_height_cm,
        points = config$arena_corner_names
      )
      tracking <- AddOFTZones(
        tracking, scale_center = 0.5, scale_periphery = 0.8,
        scale_corners = 0.4, points = config$arena_corner_names
      )
      tracking <- OFTAnalysis(
        tracking, points = "bodycentre", movement_cutoff = 5,
        integration_period = 5
      )

      code <- stringr::str_extract(inputFileName, "^[A-Za-z0-9]{4}")
      if (is.na(code)) warning("Could not extract a four-character animal code from ", inputFileName)
      novel_location <- metadata_lookup(novelLoc, code, "NovelLoc")
      metrics <- compute_socp_metrics(
        tracking, novel_location, config$fps,
        contact_distance = config$contact_distance_cm,
        proximity_range = config$proximity_range_cm,
        threshold_unit = "cm",
        min_bout_s = config$min_contact_bout_s,
        max_gap_s = config$max_contact_gap_s
      )

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
          noseObservedFraction = trackingQC$observed_fraction[trackingQC$landmark == "nose"],
          noseInterpolatedFraction = trackingQC$interpolated_fraction[trackingQC$landmark == "nose"],
          noseInvalidFraction = trackingQC$invalid_fraction[trackingQC$landmark == "nose"],
          noseLongestGapSeconds = trackingQC$longest_invalid_gap_s[trackingQC$landmark == "nose"]
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
