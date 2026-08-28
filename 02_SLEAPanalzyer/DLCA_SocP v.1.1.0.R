#' @title DLCA_SocP v1.1.0.R
#' @description Batch analysis for the social-preference workflow.
#' @version 1.1.0 (Phase 1 correctness fixes)
#'
#' Configuration is read from a YAML file, not from this script:
#'
#'   SLEAP_ANALYZER_CONFIG=my_socp.yaml Rscript "DLCA_SocP v.1.1.0.R"
#'
#' See config/socp.example.yaml for a documented template.

required_packages <- c("sp", "ggplot2", "stringr", "yaml")
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
source(file.path(script_dir, "core", "assay_config.R"))

config <- load_assay_config(resolve_config_path(assay = "SocP", script_dir = script_dir), "SocP")
repo_dir <- dirname(script_dir)

animalIDCode <- if (is.null(config$animal_id_code_file)) {
  empty_metadata_table(c("Code", "ID"))
} else {
  read_metadata_table(config$animal_id_code_file, c("Code", "ID"), if_missing = "empty")
}

batches <- if (is.null(config$batches)) "" else config$batches

for (batch in batches) {
  for (socpPhase in config$phases) {
    parts <- c(batch, socpPhase)
    parts <- parts[nzchar(parts)]
    inputDir <- do.call(file.path, as.list(c(config$input_dir, parts)))
    outputDir <- do.call(file.path, as.list(c(config$output_dir, parts)))
    metadataDir <- if (is.null(config$metadata_dir)) inputDir else config$metadata_dir
    novelLocPath <- file.path(
      metadataDir,
      paste0(config$novel_location_file_prefix, socpPhase, ".txt")
    )
    plotDir <- file.path(outputDir, "plots")
    dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

    novelLoc <- read_metadata_table(
      novelLocPath, c("Code", "NovelLoc"), sep = "\t", if_missing = "empty"
    )
    fileList <- list.files(path = inputDir, pattern = "\\.csv$", full.names = TRUE)
    dfList <- list()
    qcList <- list()

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
        tracking, points = "bodycentre",
        movement_cutoff = config$movement_cutoff_cm_s,
        integration_period = config$integration_period_frames
      )

      qcReport <- tracking_qc_report(
        tracking,
        landmarks = c("nose", "bodycentre"),
        required_landmarks = c("nose", "bodycentre"),
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
          inputFileName, " failed tracking QC: ",
          paste(qcDecision$reasons, collapse = "; "),
          ". The result is still written; exclusion is an explicit decision."
        )
      }

      code <- stringr::str_extract(inputFileName, "^[A-Za-z0-9]{4}")
      if (is.na(code)) warning("Could not extract a four-character animal code from ", inputFileName)
      novel_location <- metadata_lookup(novelLoc, code, "NovelLoc")
      metrics <- compute_socp_metrics(
        tracking, novel_location, config$fps,
        contact_distance = config$contact_distance_cm,
        body_exclusion_distance = config$body_exclusion_distance_cm,
        proximity_range = config$proximity_range_cm,
        threshold_unit = "cm",
        require_orientation = isTRUE(config$require_orientation),
        min_bout_s = config$min_bout_s,
        max_gap_s = config$max_gap_s
      )

      df <- cbind(
        data.frame(
          file = inputFileName,
          ID = metadata_lookup(animalIDCode, code, "ID"),
          Code = code,
          batch = batch,
          phase = socpPhase,
          stringsAsFactors = FALSE
        ),
        metrics$summary,
        data.frame(
          distance = tracking$Report$bodycentre.raw.distance,
          stationary = tracking$Report$bodycentre.time.stationary,
          speedMoving = tracking$Report$bodycentre.speed.moving,
          speedRaw = tracking$Report$bodycentre.raw.speed,
          qcPass = qcDecision$pass,
          qcReasons = paste(qcDecision$reasons, collapse = "; "),
          stringsAsFactors = FALSE
        ),
        qc_summary_row(qcReport)[, c(
          "valid_time_s", "valid_fraction", "observed_time_s",
          "interpolated_time_s", "longest_invalid_gap_s"
        )]
      )

      utils::write.csv(
        df, file.path(outputDir, paste0(inputFileName, "_output.csv")),
        row.names = FALSE
      )
      dfList[[length(dfList) + 1L]] <- df
      qcList[[length(qcList) + 1L]] <- qc_summary_row(qcReport)

      plots <- PlotDensityPaths(tracking, points = "bodycentre")
      if (!is.null(plots$bodycentre)) {
        ggplot2::ggsave(
          file.path(plotDir, paste0(inputFileName, "_DensityPath.png")),
          plot = plots$bodycentre, width = 7, height = 4
        )
      }
      message("Processed file ", basename(inputFile))
    }

    if (length(dfList) == 0) {
      warning("No SocP CSV files found for ", batch, " ", socpPhase, " in ", inputDir)
      next
    }
    dfCombined <- do.call(rbind, dfList)
    utils::write.csv(dfCombined, file.path(outputDir, "combined_output.csv"), row.names = FALSE)
    utils::write.csv(
      do.call(rbind, qcList), file.path(outputDir, "tracking_qc.csv"), row.names = FALSE
    )

    if (isTRUE(config$write_manifest)) {
      manifest <- run_manifest(
        config = config,
        inputs = fileList,
        packages = required_packages,
        repo_dir = repo_dir,
        extra = list(assay = "SocP", batch = batch, phase = socpPhase,
                     files = length(fileList))
      )
      assert_reproducible_run(manifest)
      write_run_manifest(manifest, file.path(outputDir, "run_manifest.yaml"))
    }
    message("Processing complete for ", batch, " ", socpPhase)
  }
}
