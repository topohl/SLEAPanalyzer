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

# Time on top of a cup is not contact (see detect_on_object()). The thresholds
# belong to the apparatus and camera, so none of them has a default.
on_cup_fields <- c("on_cup_footprint_radius_cm", "on_cup_beyond_cm", "on_cup_max_distance_cm",
                   "on_cup_min_apparent_size", "on_cup_merge_gap_s", "on_cup_min_episode_s",
                   "on_cup_size_window_s")
if (isTRUE(config$on_cup_exclusion)) {
  unset <- on_cup_fields[vapply(on_cup_fields, function(f) is.null(config[[f]]), logical(1))]
  if (length(unset) > 0) stop("on_cup_exclusion is on but unset: ", paste(unset, collapse = ", "))
}
# Landmark repair (repair_landmarks()): every setting must be stated in the configuration.
repair_fields <- c("landmark_repair_landmarks", "landmark_repair_targets", "landmark_repair_margin_cm",
                   "landmark_repair_margin_fraction", "landmark_repair_min_share", "landmark_repair_min_landmarks",
                   "landmark_repair_max_anchor_error_cm", "landmark_repair_outside_cm", "landmark_repair_pose_decay_s")
if (isTRUE(config$landmark_repair)) {
  unset <- repair_fields[vapply(repair_fields, function(f) is.null(config[[f]]), logical(1))]
  if (length(unset) > 0) stop("landmark_repair is on but unset: ", paste(unset, collapse = ", "))
  if (!all(c("nose", "bodycentre") %in% config$landmark_repair_targets)) stop("landmark_repair_targets must include nose and bodycentre")
  # a rebuilt coordinate keeps the likelihood of the point it replaced, which a cutoff would then judge
  if (!is.null(config$likelihood_cutoff)) stop("landmark_repair and likelihood_cutoff cannot be combined")
}
head_fields <- c("head_length_restore_below")
if (isTRUE(config$head_length_restore)) {
  unset <- head_fields[vapply(head_fields, function(f) is.null(config[[f]]), logical(1))]
  if (length(unset) > 0) stop("head_length_restore is on but unset: ", paste(unset, collapse = ", "))
}
if (!is.null(config$qc_decisions_file) && !file.exists(config$qc_decisions_file)) {
  stop("qc_decisions_file not found: ", config$qc_decisions_file)
}
qcDecisions <- read_qc_decisions(config$qc_decisions_file)

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
    episodeList <- list()
    repairList <- list()
    repairEpisodeList <- list()
    repairModelList <- list()

    for (inputFile in fileList) {
      inputFileName <- tools::file_path_sans_ext(basename(inputFile))
      tracking <- read_tracking_csv(file = inputFile, fps = config$fps)
      code <- stringr::str_extract(inputFileName, "^[A-Za-z0-9]{4}")

      # Calibrated first, so that the repair thresholds are in cm; interpolation
      # is linear, so its result does not depend on the order.
      tracking <- CalibrateTrackingData(
        tracking, method = "area",
        in.metric = config$arena_width_cm * config$arena_height_cm,
        points = config$arena_corner_names
      )
      repairRow <- data.frame(landmarksFlaggedS = 0, skeletonCollapseS = 0, noseReconstructedS = 0, noseRefusedS = 0,
                              bodycentreReconstructedS = 0, noseAnchorErrorP95 = NA_real_)
      if (isTRUE(config$landmark_repair)) {
        tracking <- repair_landmarks(
          tracking,
          landmarks = config$landmark_repair_landmarks,
          targets = config$landmark_repair_targets,
          margin = config$landmark_repair_margin_cm,
          margin_fraction = config$landmark_repair_margin_fraction,
          min_share = config$landmark_repair_min_share,
          min_landmarks = config$landmark_repair_min_landmarks,
          max_anchor_error = config$landmark_repair_max_anchor_error_cm,
          corners = config$arena_corner_names,
          max_outside = config$landmark_repair_outside_cm,
          pose_decay_s = config$landmark_repair_pose_decay_s
        )
        rp <- tracking$repair
        s <- rp$summary
        repairRow <- data.frame(
          landmarksFlaggedS = sum(rp$frames_flagged & !rp$collapse) / config$fps,
          skeletonCollapseS = sum(rp$collapse) / config$fps,
          noseReconstructedS = s$reconstructed_s[s$landmark == "nose"],
          noseRefusedS = s$refused_s[s$landmark == "nose"],
          bodycentreReconstructedS = s$reconstructed_s[s$landmark == "bodycentre"],
          noseAnchorErrorP95 = if (is.null(rp$models$nose)) NA_real_ else rp$models$nose$error_p95[1]
        )
        tag <- data.frame(file = inputFileName, Code = code, batch = batch, phase = socpPhase, stringsAsFactors = FALSE)
        repairList[[length(repairList) + 1L]] <- cbind(tag, s)
        if (!is.null(rp$episodes) && nrow(rp$episodes) > 0) repairEpisodeList[[length(repairEpisodeList) + 1L]] <- cbind(tag, rp$episodes)
        for (target in names(rp$models)) {
          m <- head(rp$models[[target]], 3)
          repairModelList[[length(repairModelList) + 1L]] <- cbind(tag, rank = seq_len(nrow(m)), m)
        }
      }

      # Bounded interpolation: SocP previously performed none at all, so
      # untracked frames reached the contact test as NA and were scored as
      # confident absence of social contact.
      tracking <- interpolate_tracking(
        tracking,
        landmarks = c("nose", "bodycentre"),
        max_gap_s = config$max_interpolation_gap_s,
        likelihood_cutoff = config$likelihood_cutoff
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

      onCup <- NULL
      if (isTRUE(config$on_cup_exclusion)) {
        size <- apparent_body_size(tracking)
        cups <- lapply(c(left = "socl", right = "socr"), function(cup) {
          found <- detect_on_object(
            tracking, cup, config$arena_corner_names,
            footprint_radius = config$on_cup_footprint_radius_cm,
            beyond = config$on_cup_beyond_cm,
            max_distance = config$on_cup_max_distance_cm,
            min_apparent_size = config$on_cup_min_apparent_size,
            merge_gap_s = config$on_cup_merge_gap_s,
            min_episode_s = config$on_cup_min_episode_s,
            size_window_s = config$on_cup_size_window_s,
            size = size
          )
          found$on_object <- apply_qc_decisions(
            found$on_object, qcDecisions, inputFileName, cup, config$fps
          )
          found
        })
        onCup <- list(left = cups$left$on_object, right = cups$right$on_object)
        episodes <- rbind(cups$left$episodes, cups$right$episodes)
        if (nrow(episodes) > 0) {
          episodeList[[length(episodeList) + 1L]] <- cbind(
            data.frame(file = inputFileName, Code = code, batch = batch, phase = socpPhase,
                       stringsAsFactors = FALSE),
            episodes
          )
        }
      }

      # After the on-cup detection, whose size measure uses the nose, so that
      # only the contact and proximity measures see the restored nose.
      headRow <- data.frame(headLengthCm = NA_real_, noseRestoredS = NA_real_)
      if (isTRUE(config$head_length_restore)) {
        tracking <- restore_head_length(tracking, below = config$head_length_restore_below)
        hl <- tracking$head_length
        headRow <- data.frame(headLengthCm = hl$median,
                              noseRestoredS = if (is.finite(hl$median)) sum(hl$moved) / config$fps else NA_real_)
      }

      metrics <- compute_socp_metrics(
        tracking, novel_location, config$fps,
        contact_distance = config$contact_distance_cm,
        body_exclusion_distance = config$body_exclusion_distance_cm,
        proximity_range = config$proximity_range_cm,
        threshold_unit = "cm",
        require_orientation = isTRUE(config$require_orientation),
        min_bout_s = config$min_bout_s,
        max_gap_s = config$max_gap_s,
        on_cup = onCup
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
        )],
        repairRow,
        headRow
      )

      utils::write.csv(
        df, file.path(outputDir, paste0(inputFileName, "_output.csv")),
        row.names = FALSE
      )
      dfList[[length(dfList) + 1L]] <- df
      qcList[[length(qcList) + 1L]] <- qc_summary_row(qcReport)

      # where the body centre spent its time, with each cup and its contact distance
      plots <- PlotDensityPaths(tracking, points = "bodycentre", Title = gsub("_", " ", inputFileName),
                                marks = c(socl = config$contact_distance_cm, socr = config$contact_distance_cm))
      if (!is.null(plots$bodycentre)) {
        ggplot2::ggsave(
          file.path(plotDir, paste0(inputFileName, "_DensityPath.png")),
          plot = plots$bodycentre, width = 6, height = 3.4, dpi = 300, bg = "white"
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
    if (isTRUE(config$on_cup_exclusion)) {
      # Every candidate episode, accepted or not, so each decision can be reviewed.
      episodes <- if (length(episodeList) == 0) {
        cbind(data.frame(file = character(), Code = character(), batch = character(),
                         phase = character()), empty_object_episode_table())
      } else do.call(rbind, episodeList)
      utils::write.csv(episodes, file.path(outputDir, "on_cup_episodes.csv"), row.names = FALSE)
    }
    if (isTRUE(config$landmark_repair)) {
      # Per recording and landmark, every episode, and the three best anchor models per rebuilt landmark.
      bind_or_empty <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
      utils::write.csv(bind_or_empty(repairList), file.path(outputDir, "landmark_repairs.csv"), row.names = FALSE)
      utils::write.csv(bind_or_empty(repairEpisodeList), file.path(outputDir, "landmark_repair_episodes.csv"), row.names = FALSE)
      utils::write.csv(bind_or_empty(repairModelList), file.path(outputDir, "landmark_repair_models.csv"), row.names = FALSE)
    }

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
