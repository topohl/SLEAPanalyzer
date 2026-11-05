# ============================================================================
# 28_cup_geometry.R
#
# Calibration and cup geometry of every social-preference recording, as the
# tracking sees them (no manual input): pixels per cm, the arena's aspect
# ratio (44 x 24 cm: 1.833), the distance between the cup points, the jitter
# of the cup points, the head and nose-to-tail-base length, and the apparent
# cup radius (2nd and 5th percentile of the nose distance to the cup point on
# the open side of the cup, within 60 degrees of the direction to the other
# cup, observed noses within 10 cm). Batches differ in camera setup; the cups
# are the same, so their apparent radius checks that one contact distance
# means the same in every batch.
#
#   Rscript 28_cup_geometry.R      writes results/socp_cup_geometry.csv
# ============================================================================

suppressMessages({ library(sp); library(stringr); library(yaml); library(data.table) })
SCRIPTS <- local({ a <- commandArgs(trailingOnly = FALSE); p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/") })
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
REPO <- normalizePath(file.path(PROJ, "..", ".."), winslash = "/")
Sys.setenv(SLEAP_ANALYZER_REPO_ROOT = REPO)
source(file.path(REPO, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "core", "assay_config.R"))
OUT <- file.path(PROJ, "results", "socp_cup_geometry.csv")
STAGE <- file.path(Sys.getenv("EXP9_STAGING", "C:/Users/topohl/Documents/exp9_sleap_staging"), "SocP")
config <- load_assay_config(file.path(PROJ, "config", "socp_all.yaml"), "SocP")
res <- list()
for (batch in config$batches) for (ph in c("HAB", "S1", "S2")) {
  for (f in list.files(file.path(STAGE, batch, ph), "[.]csv$", full.names = TRUE)) {
    name <- tools::file_path_sans_ext(basename(f))
    tr <- read_tracking_csv(file = f, fps = config$fps)
    corners <- sapply(c("tl", "tr", "br", "bl"), function(p) c(stats::median(tr$data[[p]]$x, na.rm = TRUE), stats::median(tr$data[[p]]$y, na.rm = TRUE)))
    side <- function(a, b) sqrt(sum((corners[, a] - corners[, b])^2))
    aspect <- (side("tl", "tr") + side("bl", "br")) / (side("tl", "bl") + side("tr", "br"))
    tr <- CalibrateTrackingData(tr, method = "area", in.metric = config$arena_width_cm * config$arena_height_cm, points = config$arena_corner_names)
    tr <- repair_landmarks(tr, landmarks = config$landmark_repair_landmarks, targets = config$landmark_repair_targets,
      margin = config$landmark_repair_margin_cm, margin_fraction = config$landmark_repair_margin_fraction, min_share = config$landmark_repair_min_share,
      min_landmarks = config$landmark_repair_min_landmarks, max_anchor_error = config$landmark_repair_max_anchor_error_cm,
      corners = config$arena_corner_names, max_outside = config$landmark_repair_outside_cm)
    nx <- tr$data$nose$x; ny <- tr$data$nose$y; ok <- as.character(tr$data$nose$status) == "observed"
    ex <- (tr$data$leftEar$x + tr$data$rightEar$x) / 2; ey <- (tr$data$leftEar$y + tr$data$rightEar$y) / 2
    hl <- stats::median(sqrt((nx - ex)^2 + (ny - ey)^2)[ok], na.rm = TRUE)
    bl <- stats::median(sqrt((nx - tr$data$tailbase$x)^2 + (ny - tr$data$tailbase$y)^2)[ok], na.rm = TRUE)
    for (cup in c("socl", "socr")) {
      cx <- stats::median(tr$data[[cup]]$x, na.rm = TRUE); cy <- stats::median(tr$data[[cup]]$y, na.rm = TRUE)
      jit <- stats::quantile(sqrt((tr$data[[cup]]$x - cx)^2 + (tr$data[[cup]]$y - cy)^2), .99, na.rm = TRUE)
      d <- sqrt((nx - cx)^2 + (ny - cy)^2); near <- ok & is.finite(d) & d < 10
      # open side of the cup (toward the other cup, +-60 degrees): where noses reach the wire
      o <- if (cup == "socl") "socr" else "socl"
      ox <- stats::median(tr$data[[o]]$x, na.rm = TRUE); oy <- stats::median(tr$data[[o]]$y, na.rm = TRUE)
      ang <- abs(((atan2(ny - cy, nx - cx) - atan2(oy - cy, ox - cx)) * 180 / pi + 180) %% 360 - 180)
      open <- near & ang < 60
      res[[length(res) + 1]] <- data.table(file = name, Code = substr(name, 1, 4), batch, phase = ph, cup,
        px_per_cm = 1 / tr$px.to.cm, aspect = aspect, head_cm = hl, nose_tail_cm = bl, cup_jitter_p99_cm = jit,
        cup_to_cup_cm = sqrt((cx - ox)^2 + (cy - oy)^2), near_s = sum(near) / 30,
        r02 = if (sum(open) > 30) stats::quantile(d[open], .02) else NA_real_, r05 = if (sum(open) > 30) stats::quantile(d[open], .05) else NA_real_,
        r50_contact = if (sum(open & d <= 6) > 30) stats::median(d[open & d <= 6]) else NA_real_)
    }
    cat(batch, ph, name, "\n")
  }
}
fwrite(rbindlist(res), OUT)
