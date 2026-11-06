# ============================================================================
# 27_figure_stills.R
#
# Stills for the figure on the tracking of the social preference test: one
# video frame per kind of misplaced nose, with the nose as tracked and as the
# analysis uses it (rebuilt by the landmark repair, or restored to the head
# length). This script finds the coordinates (in video pixels) by running the
# analysis steps of DLCA_SocP on each recording; tools/figure_stills.py draws
# them on the frames, with the cup point and the contact distance where a cup is given.
#
#   Rscript 27_figure_stills.R <out_dir>
#   python -I tools/figure_stills.py <out_dir>/stills.csv <out_dir>
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
args <- commandArgs(TRUE); if (length(args) != 1) stop("usage: Rscript 27_figure_stills.R <out_dir>")
OUT <- normalizePath(args[1], winslash = "/", mustWork = FALSE); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
STAGE <- file.path(Sys.getenv("EXP9_STAGING", "C:/Users/topohl/Documents/exp9_sleap_staging"), "SocP")
config <- load_assay_config(file.path(PROJ, "config", "socp_all.yaml"), "SocP")
VROOT <- "//mdc-berlin.net/fs/RD_Hoernberg/Tobi/E9_SIS"
VDIR <- c(B1 = "B1/20221117_SocP/cropped", B4 = "B4/SocP/cropped")
video_path <- function(batch, code, phase) {
  pat <- if (batch == "B1") sprintf("^%s_GX.*_%s[.]mp4$", code, phase) else sprintf("^E9_SIS_%s_%s_%s[.]mp4$", batch, code, phase)
  hit <- list.files(file.path(VROOT, VDIR[[batch]]), pat, full.names = TRUE); stopifnot(length(hit) == 1); hit }

# The analysis steps of DLCA_SocP up to the nose it scores: raw, repaired and interpolated, restored.
prepare <- function(batch, code, phase) {
  raw <- read_tracking_csv(file = file.path(STAGE, batch, phase, paste0(code, "_", phase, ".csv")), fps = config$fps)
  tr <- CalibrateTrackingData(raw, method = "area", in.metric = config$arena_width_cm * config$arena_height_cm, points = config$arena_corner_names)
  tr <- repair_landmarks(tr, landmarks = config$landmark_repair_landmarks, targets = config$landmark_repair_targets,
    margin = config$landmark_repair_margin_cm, margin_fraction = config$landmark_repair_margin_fraction, min_share = config$landmark_repair_min_share,
    min_landmarks = config$landmark_repair_min_landmarks, max_anchor_error = config$landmark_repair_max_anchor_error_cm,
    corners = config$arena_corner_names, max_outside = config$landmark_repair_outside_cm, pose_decay_s = config$landmark_repair_pose_decay_s)
  tr <- interpolate_tracking(tr, landmarks = c("nose", "bodycentre"), max_gap_s = config$max_interpolation_gap_s, likelihood_cutoff = config$likelihood_cutoff)
  list(raw = raw, tr = tr, rs = restore_head_length(tr, below = config$head_length_restore_below))
}

# The head under the lid: in the phase-1 recordings of batch 1 camera 2, the longest run of frames in which the restoration moves the
# nose at the lid edge to within the contact distance of a cup; its middle frame.
lid <- rbindlist(lapply(c("C7E2", "F9L3", "K4T2", "M1M3", "P9X0", "R8O3", "T2H7", "V2P0", "Y0F3"), function(code) {
  x <- prepare("B1", code, "S1"); rbindlist(lapply(c("socl", "socr"), function(cup) {
    cx <- stats::median(x$tr$data[[cup]]$x, na.rm = TRUE); cy <- stats::median(x$tr$data[[cup]]$y, na.rm = TRUE)
    d0 <- sqrt((x$tr$data$nose$x - cx)^2 + (x$tr$data$nose$y - cy)^2); d1 <- sqrt((x$rs$data$nose$x - cx)^2 + (x$rs$data$nose$y - cy)^2)
    o <- if (cup == "socl") "socr" else "socl"; ox <- stats::median(x$tr$data[[o]]$x, na.rm = TRUE)
    rx <- (x$tr$data$nose$x - cx) * sign(ox - cx); ry <- (x$tr$data$nose$y - cy) * (if (cup == "socl") 1 else -1)
    # the lid edge as in 26: beyond the contact distance, on the lid side, behind the cup point
    hit <- x$rs$head_length$moved & !is.na(d0) & d0 > config$contact_distance_cm + 0.3 & d1 <= config$contact_distance_cm & ry > 3.8 & rx < -2.5
    if (!any(hit)) return(NULL)
    r <- rle(hit); end <- cumsum(r$lengths); start <- end - r$lengths + 1; k <- which(r$values)
    data.table(Code = code, cup = cup, start = start[k], length = r$lengths[k]) })) }))
best <- lid[which.max(length)]
cat(sprintf("head under the lid: %s S1 %s, %d frames from frame %d\n", best$Code, best$cup, best$length, best$start))

# kind: what happened to the nose; the first two frames chosen from the repair review (qc/landmark_repair/review)
picks <- data.table(
  still = c("fixed_spot", "object_outside", "under_lid"),
  batch = c("B4", "B1", "B1"), Code = c("L9K1", "H9J5", best$Code), phase = c("S1", "S1", "S1"),
  frame = c(4129L, 1138L, as.integer(best$start + best$length %/% 2)), cup = c(NA, NA, best$cup),
  kind = c("nose placed on a fixed spot at the image edge", "nose placed on a tube outside the arena", "head under the cup lid, nose at the lid edge"))
rows <- list()
for (i in seq_len(nrow(picks))) {
  p <- picks[i]; x <- prepare(p$batch, p$Code, p$phase)
  f <- p$frame; k <- x$tr$px.to.cm
  if (i <= 2) {   # a rebuilt nose the head check accepted: the frame nearest to the one listed
    ok <- which(as.character(x$tr$data$nose$status) == "reconstructed"); stopifnot(length(ok) > 0)
    f <- ok[which.min(abs(ok - f))]; stopifnot(abs(f - p$frame) <= 90)
  }
  cup <- if (is.na(p$cup)) c(NA_real_, NA_real_) else c(stats::median(x$raw$data[[p$cup]]$x, na.rm = TRUE), stats::median(x$raw$data[[p$cup]]$y, na.rm = TRUE))
  tracked <- c(x$raw$data$nose$x[f], x$raw$data$nose$y[f]); used <- c(x$rs$data$nose$x[f], x$rs$data$nose$y[f]) / k
  rows[[i]] <- data.table(still = p$still, kind = p$kind, batch = p$batch, Code = p$Code, phase = p$phase, frame = f, video = video_path(p$batch, p$Code, p$phase),
    tracked_x = tracked[1], tracked_y = tracked[2], used_x = used[1], used_y = used[2], nose_status = as.character(x$rs$data$nose$status[f]),
    head_ratio = x$rs$head_length$length[f] / x$rs$head_length$median, cup_x = cup[1], cup_y = cup[2], px_per_cm = 1 / k,
    lear_x = x$tr$data$leftEar$x[f] / k, lear_y = x$tr$data$leftEar$y[f] / k, rear_x = x$tr$data$rightEar$x[f] / k, rear_y = x$tr$data$rightEar$y[f] / k,
    tracked_to_cup_cm = sqrt(sum((tracked - cup)^2)) * k, used_to_cup_cm = sqrt(sum((used - cup)^2)) * k)
}
out <- rbindlist(rows)
stopifnot("the first two noses were rebuilt by the repair" = all(out$nose_status[1:2] == "reconstructed"),
          "the third was restored into contact" = out$nose_status[3] == "reconstructed" && out$head_ratio[3] < config$head_length_restore_below &&
            out$tracked_to_cup_cm[3] > config$contact_distance_cm && out$used_to_cup_cm[3] <= config$contact_distance_cm)
fwrite(out, file.path(OUT, "stills.csv")); print(out[, .(still, Code, frame, nose_status, round(head_ratio, 2), round(tracked_to_cup_cm, 1), round(used_to_cup_cm, 1))])
