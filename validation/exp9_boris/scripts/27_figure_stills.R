# ============================================================================
# 27_figure_stills.R
#
# Stills for the figure on the tracking of the social preference test: one
# video frame per kind of misplaced nose, with the nose as tracked and as the
# analysis uses it (rebuilt by the landmark repair, or restored to the head
# length). This script finds the coordinates (in video pixels) by running the
# analysis steps of DLCA_SocP on each recording; tools/figure_stills.py draws
# them on the frames.
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

# kind: what happened to the nose; frames chosen from the repair review (qc/landmark_repair/review) and 26
picks <- data.table(
  still = c("fixed_spot", "object_outside", "under_lid"),
  batch = c("B4", "B1", "B1"), Code = c("L9K1", "H9J5", "K4T2"), phase = c("S1", "S1", "S1"), frame = c(4129L, 1138L, 2776L),
  cup = c(NA, NA, "socr"),
  kind = c("nose placed on a fixed spot at the image edge", "nose placed on a tube outside the arena", "head under the cup lid, nose at the lid edge"))
rows <- list()
for (i in seq_len(nrow(picks))) {
  p <- picks[i]
  src <- file.path(STAGE, p$batch, p$phase, paste0(p$Code, "_", p$phase, ".csv"))
  raw <- read_tracking_csv(file = src, fps = config$fps)
  tr <- CalibrateTrackingData(raw, method = "area", in.metric = config$arena_width_cm * config$arena_height_cm, points = config$arena_corner_names)
  tr <- repair_landmarks(tr, landmarks = config$landmark_repair_landmarks, targets = config$landmark_repair_targets,
    margin = config$landmark_repair_margin_cm, margin_fraction = config$landmark_repair_margin_fraction, min_share = config$landmark_repair_min_share,
    min_landmarks = config$landmark_repair_min_landmarks, max_anchor_error = config$landmark_repair_max_anchor_error_cm,
    corners = config$arena_corner_names, max_outside = config$landmark_repair_outside_cm)
  tr <- interpolate_tracking(tr, landmarks = c("nose", "bodycentre"), max_gap_s = config$max_interpolation_gap_s, likelihood_cutoff = config$likelihood_cutoff)
  rs <- restore_head_length(tr)
  f <- p$frame; k <- tr$px.to.cm
  st <- as.character(tr$data$nose$status[f])
  used <- c(x = rs$data$nose$x[f] / k, y = rs$data$nose$y[f] / k)
  rows[[i]] <- data.table(still = p$still, kind = p$kind, batch = p$batch, Code = p$Code, phase = p$phase, frame = f, video = video_path(p$batch, p$Code, p$phase),
    tracked_x = raw$data$nose$x[f], tracked_y = raw$data$nose$y[f], used_x = used[["x"]], used_y = used[["y"]], nose_status = st,
    head_ratio = rs$head_length$length[f] / rs$head_length$median,
    cup_x = if (is.na(p$cup)) NA_real_ else stats::median(raw$data[[p$cup]]$x, na.rm = TRUE),
    cup_y = if (is.na(p$cup)) NA_real_ else stats::median(raw$data[[p$cup]]$y, na.rm = TRUE), px_per_cm = 1 / k)
}
out <- rbindlist(rows)
stopifnot("the first two noses were rebuilt by the repair" = all(out$nose_status[1:2] == "reconstructed"),
          "the third was restored to the head length" = out$nose_status[3] == "observed" && out$head_ratio[3] < .75)
fwrite(out, file.path(OUT, "stills.csv")); print(out[, .(still, Code, frame, nose_status, round(head_ratio, 2))])
