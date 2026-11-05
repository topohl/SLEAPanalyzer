# ============================================================================
# 29_head_restore_impact.R
#
# What the head length restoration changes in every social-preference
# recording (all batches, habituation and phases 1 and 2): contact with each
# cup, mapped to the stranger and the other cup, without and with the
# restoration as configured, the median head length and the time restored.
# Report only; DLCA_SocP applies the restoration itself.
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
res <- list()
for (cfg in c("socp_all.yaml", "socp_hab_all.yaml")) {
  config <- load_assay_config(file.path(PROJ, "config", cfg), "SocP")
  stopifnot(isTRUE(config$head_length_restore))
  qcDecisions <- read_qc_decisions(config$qc_decisions_file)
  for (batch in config$batches) for (ph in config$phases) {
    nl <- read_metadata_table(file.path(STAGE, batch, ph, paste0(config$novel_location_file_prefix, ph, ".txt")), c("Code", "NovelLoc"), sep = "\t")
    for (f in list.files(file.path(STAGE, batch, ph), "[.]csv$", full.names = TRUE)) {
      name <- tools::file_path_sans_ext(basename(f)); code <- substr(name, 1, 4)
      tr <- read_tracking_csv(file = f, fps = config$fps)
      tr <- CalibrateTrackingData(tr, method = "area", in.metric = config$arena_width_cm * config$arena_height_cm, points = config$arena_corner_names)
      tr <- repair_landmarks(tr, landmarks = config$landmark_repair_landmarks, targets = config$landmark_repair_targets,
        margin = config$landmark_repair_margin_cm, margin_fraction = config$landmark_repair_margin_fraction, min_share = config$landmark_repair_min_share,
        min_landmarks = config$landmark_repair_min_landmarks, max_anchor_error = config$landmark_repair_max_anchor_error_cm,
        corners = config$arena_corner_names, max_outside = config$landmark_repair_outside_cm)
      tr <- interpolate_tracking(tr, landmarks = c("nose", "bodycentre"), max_gap_s = config$max_interpolation_gap_s, likelihood_cutoff = config$likelihood_cutoff)
      size <- apparent_body_size(tr)
      on <- lapply(c(left = "socl", right = "socr"), function(cup) {
        x <- detect_on_object(tr, cup, config$arena_corner_names, footprint_radius = config$on_cup_footprint_radius_cm, beyond = config$on_cup_beyond_cm,
          max_distance = config$on_cup_max_distance_cm, min_apparent_size = config$on_cup_min_apparent_size, merge_gap_s = config$on_cup_merge_gap_s,
          min_episode_s = config$on_cup_min_episode_s, size_window_s = config$on_cup_size_window_s, size = size)
        apply_qc_decisions(x$on_object, qcDecisions, name, cup, config$fps) })
      loc <- metadata_lookup(nl, code, "NovelLoc")
      sc <- function(t) compute_socp_metrics(t, loc, config$fps, contact_distance = config$contact_distance_cm,
        body_exclusion_distance = config$body_exclusion_distance_cm, proximity_range = config$proximity_range_cm, threshold_unit = "cm", on_cup = on)$summary
      rs <- restore_head_length(tr, below = config$head_length_restore_below)
      b <- sc(tr); a <- sc(rs)
      res[[length(res) + 1]] <- data.table(file = name, Code = code, batch, phase = ph, novelLoc = loc,
        left_before = b$contactLeft, right_before = b$contactRight, left_after = a$contactLeft, right_after = a$contactRight,
        stranger_before = b$contactNovel, other_before = b$contactFamiliar, stranger_after = a$contactNovel, other_after = a$contactFamiliar,
        head_length_cm = rs$head_length$median, nose_restored_s = sum(rs$head_length$moved) / config$fps, total_s = length(rs$head_length$moved) / config$fps)
    }
  }
}
out <- rbindlist(res)
stopifnot(nrow(out) == 348, !anyNA(out$novelLoc))
fwrite(out, OUT)
cat(nrow(out), "recordings; contact gain (s) per cup: median", round(stats::median(c(out$left_after - out$left_before, out$right_after - out$right_before)), 2), "\n")
