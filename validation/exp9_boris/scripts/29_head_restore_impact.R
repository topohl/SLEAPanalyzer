# ============================================================================
# 29_head_restore_impact.R
#
# What the landmark repair as now configured (head check, pose of the moment)
# and the head length restoration change in every social-preference recording
# (all batches, habituation and phases 1 and 2): contact with each cup, mapped
# to the stranger and the other cup, (1) as released (repair without the head
# check and the pose, no restoration), (2) repaired as configured, (3) repaired
# and restored; the nose time rebuilt and refused, the time restored, the
# median head length and the longest stretch without nose or body centre.
# Report only; DLCA_SocP applies all of it itself.
#
#   Rscript 29_head_restore_impact.R      writes results/socp_head_restore_impact.csv
# ============================================================================

suppressMessages({ library(sp); library(stringr); library(yaml); library(data.table) })
SCRIPTS <- local({ a <- commandArgs(trailingOnly = FALSE); p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/") })
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
REPO <- normalizePath(file.path(PROJ, "..", ".."), winslash = "/")
Sys.setenv(SLEAP_ANALYZER_REPO_ROOT = REPO)
source(file.path(REPO, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "Behavioral_Metrics_Phase1.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "core", "assay_config.R"))
OUT <- file.path(PROJ, "results", "socp_head_restore_impact.csv")
STAGE <- file.path(Sys.getenv("EXP9_STAGING", "C:/Users/topohl/Documents/exp9_sleap_staging"), "SocP")
longest <- function(x) { r <- rle(x); if (any(r$values)) max(r$lengths[r$values]) else 0L }
res <- list()
for (cfg in c("socp_all.yaml", "socp_hab_all.yaml")) {
  config <- load_assay_config(file.path(PROJ, "config", cfg), "SocP")
  stopifnot(isTRUE(config$head_length_restore))
  qcDecisions <- read_qc_decisions(config$qc_decisions_file)
  for (batch in config$batches) for (ph in config$phases) {
    nl <- read_metadata_table(file.path(STAGE, batch, ph, paste0(config$novel_location_file_prefix, ph, ".txt")), c("Code", "NovelLoc"), sep = "\t")
    for (f in list.files(file.path(STAGE, batch, ph), "[.]csv$", full.names = TRUE)) {
      name <- tools::file_path_sans_ext(basename(f)); code <- substr(name, 1, 4); loc <- metadata_lookup(nl, code, "NovelLoc")
      raw <- read_tracking_csv(file = f, fps = config$fps)
      raw <- CalibrateTrackingData(raw, method = "area", in.metric = config$arena_width_cm * config$arena_height_cm, points = config$arena_corner_names)
      run <- function(as_released, restore) {
        tr <- repair_landmarks(raw, landmarks = config$landmark_repair_landmarks, targets = config$landmark_repair_targets,
          margin = config$landmark_repair_margin_cm, margin_fraction = config$landmark_repair_margin_fraction, min_share = config$landmark_repair_min_share,
          min_landmarks = config$landmark_repair_min_landmarks, max_anchor_error = config$landmark_repair_max_anchor_error_cm,
          corners = config$arena_corner_names, max_outside = config$landmark_repair_outside_cm, head_check = !as_released,
          pose_decay_s = if (as_released) 0 else config$landmark_repair_pose_decay_s)
        s <- tr$repair$summary
        tr <- interpolate_tracking(tr, landmarks = c("nose", "bodycentre"), max_gap_s = config$max_interpolation_gap_s, likelihood_cutoff = config$likelihood_cutoff)
        size <- apparent_body_size(tr)
        on <- lapply(c(left = "socl", right = "socr"), function(cup) {
          x <- detect_on_object(tr, cup, config$arena_corner_names, footprint_radius = config$on_cup_footprint_radius_cm, beyond = config$on_cup_beyond_cm,
            max_distance = config$on_cup_max_distance_cm, min_apparent_size = config$on_cup_min_apparent_size, merge_gap_s = config$on_cup_merge_gap_s,
            min_episode_s = config$on_cup_min_episode_s, size_window_s = config$on_cup_size_window_s, size = size)
          apply_qc_decisions(x$on_object, qcDecisions, name, cup, config$fps) })
        restored <- NA_real_; head_cm <- NA_real_
        if (restore) { tr <- restore_head_length(tr, below = config$head_length_restore_below)
          restored <- sum(tr$head_length$moved) / config$fps; head_cm <- tr$head_length$median }
        m <- compute_socp_metrics(tr, loc, config$fps, contact_distance = config$contact_distance_cm,
          body_exclusion_distance = config$body_exclusion_distance_cm, proximity_range = config$proximity_range_cm, threshold_unit = "cm", on_cup = on)$summary
        gap <- !(is.finite(tr$data$nose$x) & is.finite(tr$data$bodycentre$x))
        list(m = m, rebuilt = s$reconstructed_s[s$landmark == "nose"], refused = s$refused_s[s$landmark == "nose"], restored = restored, head_cm = head_cm,
             gap_s = longest(gap) / config$fps, n = length(gap))
      }
      a <- run(TRUE, FALSE); b <- run(FALSE, FALSE); c3 <- run(FALSE, TRUE)
      res[[length(res) + 1]] <- data.table(file = name, Code = code, batch, phase = ph, novelLoc = loc,
        left_released = a$m$contactLeft, right_released = a$m$contactRight, stranger_released = a$m$contactNovel, other_released = a$m$contactFamiliar,
        left_repaired = b$m$contactLeft, right_repaired = b$m$contactRight, stranger_repaired = b$m$contactNovel, other_repaired = b$m$contactFamiliar,
        left_final = c3$m$contactLeft, right_final = c3$m$contactRight, stranger_final = c3$m$contactNovel, other_final = c3$m$contactFamiliar,
        nose_rebuilt_released_s = a$rebuilt, nose_rebuilt_s = b$rebuilt, nose_refused_s = b$refused, nose_restored_s = c3$restored, head_length_cm = c3$head_cm,
        longest_gap_released_s = a$gap_s, longest_gap_s = b$gap_s, total_s = c3$n / config$fps)
    }
  }
}
out <- rbindlist(res)
stopifnot(nrow(out) == 348, !anyNA(out$novelLoc))
fwrite(out, OUT)
cat(nrow(out), "recordings; nose rebuilt", round(sum(out$nose_rebuilt_s)), "s, refused by the head check", round(sum(out$nose_refused_s)), "s; longest gap above 5 s in",
    sum(out$longest_gap_s > 5), "recordings (released:", sum(out$longest_gap_released_s > 5), ")\n")
