# ============================================================================
# 25_repair_review.R
#
# Review material for the landmark repair of the Exp9 social preference test:
# picks episodes from the SocP outputs, writes each picked recording's
# tracking after the repair (same layout as the staged file, nose and body
# centre replaced, before interpolation) and the episode and video lists for
# tools/qc_review.py, which then draws the original nose and body centre (red,
# green) and the repaired ones (cyan, magenta rings) on stills and clips.
#
#   Rscript 25_repair_review.R <socp_output_root> <work_dir>
#   python -I tools/qc_review.py <work_dir>/episodes.csv <work_dir>/videos.csv <out_dir>
#
# Picked: the longest rebuilt-nose episodes (at least 1 s, mostly rebuilt) of
# the 12 recordings with the most rebuilt nose, three each; the longest of
# each named case (H9J5 S1 nose on a tube, Q5M7 reflection, C1U2 fixed spot);
# every collapsed-skeleton episode of at least 1 s; the ten longest split
# skeletons. work_dir holds the repaired tracking files (about 10 MB each).
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
args <- commandArgs(TRUE); if (length(args) != 2) stop("usage: Rscript 25_repair_review.R <socp_output_root> <work_dir>")
ROOT <- args[1]; WORK <- normalizePath(args[2], winslash = "/", mustWork = FALSE)
dir.create(file.path(WORK, "tracking"), recursive = TRUE, showWarnings = FALSE)
STAGE <- Sys.getenv("EXP9_STAGING", "C:/Users/topohl/Documents/exp9_sleap_staging")
SOCP <- load_assay_config(file.path(PROJ, "config", "socp_all.yaml"), "SocP")
VROOT <- "//mdc-berlin.net/fs/RD_Hoernberg/Tobi/E9_SIS"
VDIR <- c(B1 = "B1/20221117_SocP/cropped", B2 = "B2/20230218_SocP/cropped", B3 = "B3/SocP/cropped", B4 = "B4/SocP/cropped",
          B5 = "B5/20231210_SocP/cropped", B6 = "B6/SocP/cropped")
ALIAS <- c("B4 B7C9 S2" = "B7C2", "B2 Q6Q6 HAB" = "O6O6", "B2 Q6Q6 S1" = "O6O6", "B2 Q6Q6 S2" = "O6O6")   # video names with a typo
video_path <- function(batch, code, phase) {
  v <- ALIAS[paste(batch, code, phase)]; if (!is.na(v)) code <- unname(v)
  pat <- if (batch == "B1") sprintf("^%s_GX.*_%s[.]mp4$", code, phase) else sprintf("^E9_SIS_%s_%s_%s[.]mp4$", batch, code, phase)
  hit <- list.files(file.path(VROOT, VDIR[[batch]]), pat, full.names = TRUE)
  if (length(hit) != 1) stop(batch, " ", code, " ", phase, ": ", length(hit), " videos")
  hit }

rd <- function(name) rbindlist(lapply(list.files(ROOT, paste0("^", name, "$"), recursive = TRUE, full.names = TRUE), function(f)
  if (file.size(f) > 5) fread(f)), fill = TRUE)
ep <- rd("landmark_repair_episodes.csv"); co <- rd("combined_output.csv")
stopifnot("repair episodes found" = nrow(ep) > 0,
          "every SocP output comes from a run with the landmark repair" = "noseReconstructedS" %in% names(co) && !anyNA(co$noseReconstructedS))
cat("SocP outputs:", nrow(co), "recordings in", uniqueN(co[, paste(batch, phase)]), "batch/phase folders\n")
nose <- ep[target == "nose" & check %in% c("inconsistent_landmark", "outside_arena") & duration_s >= 1 & reconstructed_share >= 0.5]
top <- co[order(-noseReconstructedS)][1:12, .(batch, phase, Code)]
picks <- rbind(
  nose[top, on = c("batch", "phase", "Code"), nomatch = 0][order(-duration_s), head(.SD, 3), by = .(batch, phase, Code)][, group := "Nose rebuilt: recordings with the most"],
  nose[paste(Code, phase) %in% c("H9J5 S1", "Q5M7 HAB", "Q5M7 S1", "Q5M7 S2", "C1U2 S1", "C1U2 S2")][order(-duration_s), head(.SD, 1), by = .(batch, phase, Code)][, group := "Nose rebuilt: named cases"],
  ep[check == "skeleton_collapse" & duration_s >= 1][, group := "Skeleton collapsed (left untracked)"],
  head(ep[check == "skeleton_split"][order(-duration_s)], 10)[, group := "Skeleton split (the consistent majority kept)"], fill = TRUE)
picks <- unique(picks, by = c("batch", "phase", "Code", "check", "target", "start_frame"))
picks[, episode_id := sprintf("rp_%s_%s_%s_%02d", batch, Code, phase, seq_len(.N)), by = .(batch, phase, Code)]
cat(nrow(picks), "episodes from", uniqueN(picks[, paste(batch, phase, Code)]), "recordings\n")

# the repaired tracking of each picked recording, in the staged layout (pixels), nose and body centre replaced
for (k in unique(picks[, paste(batch, phase, Code)])) {
  b <- strsplit(k, " ")[[1]]; src <- file.path(STAGE, "SocP", b[1], b[2], paste0(b[3], "_", b[2], ".csv"))
  tr <- read_tracking_csv(file = src, fps = SOCP$fps)
  tr <- CalibrateTrackingData(tr, method = "area", in.metric = SOCP$arena_width_cm * SOCP$arena_height_cm, points = SOCP$arena_corner_names)
  rp <- repair_landmarks(tr, landmarks = SOCP$landmark_repair_landmarks, targets = SOCP$landmark_repair_targets, margin = SOCP$landmark_repair_margin_cm,
                         margin_fraction = SOCP$landmark_repair_margin_fraction, min_share = SOCP$landmark_repair_min_share,
                         min_landmarks = SOCP$landmark_repair_min_landmarks, max_anchor_error = SOCP$landmark_repair_max_anchor_error_cm,
                         corners = SOCP$arena_corner_names, max_outside = SOCP$landmark_repair_outside_cm)
  head3 <- readLines(src, n = 3); body <- fread(src, skip = 3, header = FALSE)
  parts <- strsplit(head3[2], ",")[[1]]; coords <- strsplit(head3[3], ",")[[1]]
  stopifnot(nrow(body) == length(rp$data$nose$x))
  for (p in c("nose", "bodycentre")) for (xy in c("x", "y")) {
    j <- which(parts == p & coords == xy); stopifnot(length(j) == 1)
    set(body, j = j, value = rp$data[[p]][[xy]] / tr$px.to.cm) }
  out <- file.path(WORK, "tracking", paste0(b[3], "_", b[2], "_repaired.csv"))
  writeLines(head3, out); fwrite(body, out, append = TRUE, col.names = FALSE, na = "NaN")
  picks[paste(batch, phase, Code) == k, `:=`(tracking_file = src, repaired_file = out)]
}
picks[, label := fifelse(check == "skeleton_collapse", "COLLAPSED: untracked", fifelse(check == "skeleton_split", "SPLIT: majority kept",
                         sprintf("%s REBUILT (%s)", toupper(target), check)))]
ids <- unique(co[, .(Code, ID)])
picks <- ids[picks, on = "Code"]
fwrite(picks[, .(episode_id, group, ID, Code, batch, phase, check, landmark = target, start_frame, end_frame, start_s = round(start_s, 2),
                 duration_s = round(duration_s, 2), reconstructed_share = round(reconstructed_share, 2), label, tracking_file, repaired_file)],
       file.path(WORK, "episodes.csv"))
vids <- unique(picks[, .(batch, phase, Code, tracking_file)])
vids[, video_path := mapply(video_path, batch, Code, phase)]
fwrite(vids[, .(tracking_file, video_path)], file.path(WORK, "videos.csv"))
cat("wrote", file.path(WORK, c("episodes.csv", "videos.csv")), "and", nrow(vids), "repaired tracking files\n")
