# ============================================================================
# 23_plausibility_screen.R
#
# Frame-level plausibility screen of every staged Exp9 recording (social
# preference, novel object, open field, plus maze; tests and habituation
# trials). Report only: nothing here changes a score. Each recording is read,
# interpolated and calibrated as its assay script does, then the checks of
# core/plausibility.R run on it:
#
#   jump           body centre faster than max_plausible_speed_cm_s
#   skeleton       nose-tailbase length outside +/- 50% of the recording median
#   outside_arena  body centre more than 1 cm outside the arena corners
#                  (not the plus maze, whose corners are arm ends)
#   raised         apparent body size >= 1.15 (0.5 s running median) for >= 1 s:
#                  the mouse up on something, closer to the camera
#   on_object      detect_on_object() at the cups (SocP) and objects (NOR),
#                  with the contact the assay's own rule scores inside each
#                  accepted episode; SocP uses the release rule (on_cup_*),
#                  NOR the setting in NOR_ON_OBJECT below (not applied by the
#                  NOR scoring: report only)
#   repair         repair_landmarks() with the SocP settings (the plus maze
#                  without the outside-arena step): seconds with a removed
#                  landmark, with a split or collapsed skeleton, with a rebuilt
#                  nose or body centre; for the novel-object test also contact
#                  with each object before and after (report only: only SocP
#                  applies the repair)
#
# Writes results/plausibility_episodes.csv (one row per episode) and
# results/plausibility_summary.csv (one row per recording: seconds and
# episodes per check). tools/qc_review.py renders stills and clips for any
# selection of the episodes. See docs/assay_definitions.md, "Plausibility
# checks".
#
# Inputs: the staged recordings (C:/Users/topohl/Documents/exp9_sleap_staging)
# and the release configs in ../config. Runtime about an hour.
# ============================================================================

suppressMessages({ library(sp); library(ggplot2); library(stringr); library(yaml); library(data.table) })
SCRIPTS <- local({ a <- commandArgs(trailingOnly = FALSE); p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/") })
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
REPO <- normalizePath(file.path(PROJ, "..", ".."), winslash = "/")
Sys.setenv(SLEAP_ANALYZER_REPO_ROOT = REPO)
source(file.path(REPO, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "Behavioral_Metrics_Phase1.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "core", "assay_config.R"))
RES <- file.path(PROJ, "results")
STAGE <- Sys.getenv("EXP9_STAGING", "C:/Users/topohl/Documents/exp9_sleap_staging")
cfg <- function(f, assay) load_assay_config(file.path(PROJ, "config", f), assay)
SOCP <- cfg("socp_all.yaml", "SocP"); NOR <- cfg("nor_all.yaml", "NOR"); NORH <- cfg("nor_hab_all.yaml", "NOR"); EPM <- cfg("epm_all.yaml", "EPM")
CORNERS <- c("tl", "tr", "br", "bl"); FPS <- 30; MAX_SPEED <- 100

# NOR objects, report only: checked on batch 1 against the 28 manual climbing bouts (2026-11-05). The NOR scoring already
# leaves out most full climbs through its head-angle rule; this setting removes the most contact inside climbing at 96%
# precision (3.5 cm, 1.04: 16 of 28 bouts, the low rectangular object raises a mouse only ~1.02). See docs/assay_definitions.md.
NOR_ON_OBJECT <- list(footprint_radius = 3.5, beyond = NULL, max_distance = 8, min_apparent_size = 1.04, merge_gap_s = 0.5,
                      min_episode_s = 0.2, size_window_s = 0.5)

recordings <- rbindlist(list(
  data.table(assay = "SocP", phase = rep(c("HAB", "S1", "S2"), each = 6), batch = rep(paste0("B", 1:6), 3))[
    , .(file = list.files(file.path(STAGE, "SocP", batch, phase), "[.]csv$", full.names = TRUE)), by = .(assay, phase, batch)],
  data.table(assay = "NOR", phase = "NOV", batch = paste0("B", 1:6))[
    , .(file = list.files(file.path(STAGE, "NOR", batch), "[.]csv$", full.names = TRUE)), by = .(assay, phase, batch)],
  data.table(assay = "NOR", phase = "HAB", batch = paste0("B", 2:6))[
    , .(file = list.files(file.path(STAGE, "NOR_HAB", batch), "[.]csv$", full.names = TRUE)), by = .(assay, phase, batch)],
  data.table(assay = "OFT", phase = "OFT", batch = paste0("B", 1:6))[
    , .(file = list.files(file.path(STAGE, "OFT", batch), "[.]csv$", full.names = TRUE)), by = .(assay, phase, batch)],
  data.table(assay = "EPM", phase = "EPM", batch = paste0("B", 1:6))[
    , .(file = list.files(file.path(STAGE, "EPM", batch), "[.]csv$", full.names = TRUE)), by = .(assay, phase, batch)]))
recordings[, Code := str_extract(basename(file), "^[A-Za-z0-9]{4}")]
cat(nrow(recordings), "recordings\n")

prepare <- function(r, interpolate = TRUE) {
  tr <- read_tracking_csv(file = r$file, fps = FPS)
  if (interpolate) tr <- interpolate_tracking(tr, landmarks = c("nose", "bodycentre"), max_gap_s = 0.2, likelihood_cutoff = NULL)
  if (r$assay == "EPM") {
    tr <- CalibrateTrackingData(tr, method = "distance", in.metric = EPM$calibration_distance_cm, points = EPM$calibration_points)
  } else {
    side <- switch(r$assay, SocP = c(SOCP$arena_width_cm, SOCP$arena_height_cm), NOR = c(NOR$arena_width_cm, NOR$arena_height_cm), OFT = c(49, 49))
    tr <- CalibrateTrackingData(tr, method = "area", in.metric = prod(side), points = CORNERS)
  }
  if (r$assay == "NOR") tr <- hold_static_landmarks(tr, c("objL", "objR"), NOR$static_objects_min_separation_cm)
  tr
}
novel_side <- function(r) {
  p <- if (r$assay == "SocP") file.path(STAGE, "SocP", r$batch, r$phase, paste0("novelLoc", r$phase, ".txt")) else file.path(NOR$metadata_dir, r$batch, NOR$novel_location_file)
  x <- read_metadata_table(p, c("Code", "NovelLoc"), sep = "\t", if_missing = "empty"); metadata_lookup(x, r$Code, "NovelLoc")
}
episode_table <- function(r, rows) if (nrow(rows)) cbind(data.table(assay = r$assay, phase = r$phase, batch = r$batch, Code = r$Code), rows) else NULL

screen_one <- function(r) {
  tr <- prepare(r); n <- length(get_tracking_frames(tr))
  size <- apparent_body_size(tr)
  raised <- rolling_median(size, 15L) >= 1.15
  out <- list(
    jump_episodes(tr, "bodycentre", max_speed = MAX_SPEED),
    skeleton_episodes(tr, "nose", "tailbase", tolerance = 0.5),
    if (r$assay != "EPM") outside_arena_episodes(tr, "bodycentre", CORNERS, tolerance = 1),
    { e <- episode_rows(raised %in% TRUE, FPS, 0.2, check = "raised", target = "body", statistic = size,
                        statistic_name = "median_apparent_size", summary = stats::median); e[e$duration_s >= 1, , drop = FALSE] })
  out <- Filter(function(e) !is.null(e) && nrow(e) > 0, out)
  objects <- switch(r$assay, SocP = c(left = "socl", right = "socr"), NOR = c(left = "objL", right = "objR"), NULL)
  p <- if (r$assay == "SocP") list(footprint_radius = SOCP$on_cup_footprint_radius_cm, beyond = SOCP$on_cup_beyond_cm,
                                   max_distance = SOCP$on_cup_max_distance_cm, min_apparent_size = SOCP$on_cup_min_apparent_size,
                                   merge_gap_s = SOCP$on_cup_merge_gap_s, min_episode_s = SOCP$on_cup_min_episode_s,
                                   size_window_s = SOCP$on_cup_size_window_s) else NOR_ON_OBJECT
  if (!is.null(objects) && !is.null(p)) {
    loc <- novel_side(r)
    events <- if (r$assay == "SocP") compute_socp_metrics(tr, loc, FPS, contact_distance = SOCP$contact_distance_cm,
                  body_exclusion_distance = SOCP$body_exclusion_distance_cm, proximity_range = SOCP$proximity_range_cm, threshold_unit = "cm")$events
              else { Q <- if (r$phase == "HAB") NORH else NOR; compute_nor_metrics(tr, loc, FPS, contact_geometry = Q$contact_geometry, contact_distance = NOR$contact_distance_cm,
                  body_exclusion_distance = NOR$body_exclusion_distance_cm, object_box_width = NOR$object_box_width_cm,
                  object_box_height = NOR$object_box_height_cm, contact_angle = NOR$contact_angle_deg, proximity_range = NOR$proximity_range_cm,
                  proximity_angle = NOR$proximity_angle_deg)$events }
    for (side in names(objects)) {
      d <- do.call(detect_on_object, c(list(tracking = tr, object = objects[[side]], corners = CORNERS, size = size), p))$episodes
      if (!nrow(d)) next
      contact <- events[[paste0("contact_", side)]] %in% TRUE
      out[[length(out) + 1]] <- data.frame(check = "on_object", target = d$object, start_frame = d$start_frame, end_frame = d$end_frame,
        start_s = d$start_s, duration_s = d$duration_s, statistic_name = "median_apparent_size", statistic = d$median_apparent_size,
        decision = d$decision, contact_s = vapply(seq_len(nrow(d)), function(i) sum(contact[d$start_frame[i]:d$end_frame[i]]) / FPS, 0))
    }
  }
  e <- rbindlist(out, fill = TRUE)
  for (cn in c("decision", "contact_s")) if (nrow(e) && !cn %in% names(e)) e[, (cn) := NA]
  # Landmark repair with the SocP settings (report only here): what it would flag and rebuild, and for the novel-object
  # test what it would change in contact and D2. In the order of DLCA_SocP: calibrate, repair, then interpolate.
  raw <- prepare(r, interpolate = FALSE)
  rp <- tryCatch(repair_landmarks(raw, landmarks = intersect(SOCP$landmark_repair_landmarks, names(tr$data)), targets = SOCP$landmark_repair_targets,
                                  margin = SOCP$landmark_repair_margin_cm, margin_fraction = SOCP$landmark_repair_margin_fraction,
                                  min_share = SOCP$landmark_repair_min_share, min_landmarks = SOCP$landmark_repair_min_landmarks,
                                  max_anchor_error = SOCP$landmark_repair_max_anchor_error_cm,
                                  corners = if (r$assay != "EPM") CORNERS, max_outside = SOCP$landmark_repair_outside_cm), error = function(e) NULL)
  if (is.null(rp)) warning("landmark repair failed for ", basename(r$file))
  repair <- if (is.null(rp)) NULL else { rs <- rp$repair$summary
    list(repair_flagged_s = sum(rp$repair$frames_flagged & !rp$repair$collapse) / FPS, repair_split_s = sum(rp$repair$split) / FPS,
         repair_collapse_s = sum(rp$repair$collapse) / FPS,
         nose_rebuilt_s = rs$reconstructed_s[rs$landmark == "nose"], body_rebuilt_s = rs$reconstructed_s[rs$landmark == "bodycentre"]) }
  if (!is.null(rp) && r$assay == "NOR" && r$phase == "NOV") {
    sc <- function(t) { t <- interpolate_tracking(t, landmarks = c("nose", "bodycentre"), max_gap_s = NOR$max_interpolation_gap_s)
      compute_nor_metrics(t, novel_side(r), FPS, contact_geometry = NOR$contact_geometry, contact_distance = NOR$contact_distance_cm,
        body_exclusion_distance = NOR$body_exclusion_distance_cm, object_box_width = NOR$object_box_width_cm, object_box_height = NOR$object_box_height_cm,
        contact_angle = NOR$contact_angle_deg, proximity_range = NOR$proximity_range_cm, proximity_angle = NOR$proximity_angle_deg)$summary }
    b <- sc(raw); a <- sc(rp)
    repair <- c(repair, list(nor_novel_before = b$contactNov, nor_novel_after = a$contactNov, nor_familiar_before = b$contactFam, nor_familiar_after = a$contactFam))
  }
  list(episodes = episode_table(r, e), frames = n, repair = repair)
}

eps <- list(); summ <- list()
for (i in seq_len(nrow(recordings))) {
  r <- recordings[i]
  res <- tryCatch(screen_one(r), error = function(e) { message("FAILED ", r$file, ": ", conditionMessage(e)); NULL })
  if (is.null(res)) { summ[[i]] <- data.table(r[, .(assay, phase, batch, Code)], frames = NA_integer_, failed = TRUE); next }
  e <- res$episodes; eps[[i]] <- e
  s <- data.table(r[, .(assay, phase, batch, Code)], frames = res$frames, failed = FALSE)
  for (ck in c("jump", "skeleton", "outside_arena", "raised")) {
    x <- if (is.null(e)) NULL else e[check == ck]
    s[, (paste0(ck, "_n")) := if (is.null(x)) 0L else nrow(x)][, (paste0(ck, "_s")) := if (is.null(x)) 0 else round(sum(x$duration_s), 2)]
  }
  x <- if (is.null(e)) NULL else e[check == "on_object" & decision == "on_object"]
  s[, on_object_n := if (is.null(x)) 0L else nrow(x)][, on_object_s := if (is.null(x)) 0 else round(sum(x$duration_s), 2)]
  s[, on_object_contact_s := if (is.null(x)) 0 else round(sum(x$contact_s), 2)]
  if (!is.null(res$repair)) for (nm in names(res$repair)) s[, (nm) := round(res$repair[[nm]], 3)]
  summ[[i]] <- s
  if (i %% 50 == 0) cat(i, "of", nrow(recordings), "\n")
}
E <- rbindlist(eps, fill = TRUE); S <- rbindlist(summ, fill = TRUE)
fwrite(E, file.path(RES, "plausibility_episodes.csv")); fwrite(S, file.path(RES, "plausibility_summary.csv"))
cat("\nrecordings per assay and phase with at least one episode:\n")
print(S[, .(recordings = .N, failed = sum(failed), jump = sum(jump_n > 0, na.rm = TRUE), skeleton = sum(skeleton_n > 0, na.rm = TRUE),
            outside_arena = sum(outside_arena_n > 0, na.rm = TRUE), raised = sum(raised_n > 0, na.rm = TRUE),
            on_object = sum(on_object_n > 0, na.rm = TRUE), on_object_contact_s = sum(on_object_contact_s, na.rm = TRUE)), by = .(assay, phase)])
cat("wrote", file.path(RES, c("plausibility_episodes.csv", "plausibility_summary.csv")), sep = "\n  ")
