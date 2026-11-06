# ============================================================================
# 31_light_checks.R
#
# A light check of every staged Exp9 recording: is it the right trial and the
# whole trial, is the arena where it should be, does the tracked mouse move
# like a mouse, and do each test's own numbers add up. Written for the novel
# object, open field and plus maze tests, which have not had the frame-by-frame
# review the social preference test had; the social preference trials are
# included so that one table covers every test. Report only: nothing here
# changes a score. The frame-level checks of 23_plausibility_screen.R are
# joined in as information, not rerun.
#
#   length     recorded duration (frames at the cameras' 29.97 per second)
#              against the protocol: 5 min plus maze and novel object trials,
#              10 min open field and social preference trials
#   last_bin   open field: a last 60 s bin with under 90% of its frames, which
#              DLCA_OFT counted as a full minute in its habituation slope
#              before 2026-11-06
#   duplicate  the same trial twice: identical files, or two recordings of a
#              batch whose body-centre tracks coincide (median distance under
#              0.15 cm) at some offset within 60 s, e.g. a phase cut from the
#              wrong part of a video
#   landmarks  a landmark the test uses that is never, or rarely, tracked
#   arena      side ratio of the tracked arena against the real one (plus
#              maze: the width of an arm read back through the 60 cm
#              calibration); corner drift (95th percentile distance of each
#              corner from its median), the largest share of frames with a
#              corner more than 2 cm off (a corner lost while the mouse covers
#              it), and camera shift (all corners moving together, largest
#              10 s median). Zones, the cm scale and fixed objects come from
#              median positions, which only a shift or a loss in most frames
#              would move.
#   mouse      median nose-tailbase length while the mouse moves (body centre
#              faster than 10 cm/s) against its test, sex and batch, and the
#              share of moving frames in which the body centre moves toward the
#              nose (a swapped nose and tail base would reverse it)
#   spikes     the whole skeleton (nose, body centre, tail base) more than 1 cm
#              off the midpoint of the frames either side: a pose taken from
#              another video frame. Reported as the path it adds, in % of the
#              distance; body-centre-only flicks are listed for information
#   objects    novel object trials: the share of frames in which the tracked
#              objects stand apart (which sets their fixed positions) and, per
#              object, in which the tracked point is within 2 cm of its fixed
#              position; the distance between them is information only (hand
#              placement varies)
#   off_zone   plus maze: share of valid time in no zone; open field: body-
#              centre time the scoring dropped (outside the arena or a jump;
#              median 2.4 s, flagged above 15 s, 2.5% of the test)
#   contact    novel object test: total object contact, which D2 rests on
#   readouts   robust z (median and MAD within test and sex) of the main
#              readouts and distance, on the trial they come from; |z| > 3.5
#              is listed for a look, never a reason to exclude by itself. Also
#              a phase-1 social preference below 0 (a swapped stimulus side?)
#   video      tools/video_qc flags, read from results/video_qc/video_qc.csv
#              when it is there (manifest from 32_video_qc_manifest.R)
#
# Writes results/light_checks_recordings.csv (one row per recording, values
# and flags), results/light_checks_summary.csv (one row per test and trial)
# and results/light_checks_duplicates.csv (pairs of recordings that are the
# same trial). Inputs: the staged recordings, the release named in RELEASE,
# and results/plausibility_summary.csv. Runtime about 20 minutes.
# ============================================================================

suppressMessages({ library(sp); library(stringr); library(yaml); library(data.table) })
SCRIPTS <- local({ a <- commandArgs(trailingOnly = FALSE); p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/") })
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
REPO <- normalizePath(file.path(PROJ, "..", ".."), winslash = "/")
Sys.setenv(SLEAP_ANALYZER_REPO_ROOT = REPO)
source(file.path(REPO, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "core", "assay_config.R"))
RES <- file.path(PROJ, "results")
STAGE <- Sys.getenv("EXP9_STAGING", "C:/Users/topohl/Documents/exp9_sleap_staging")
RELEASE <- Sys.getenv("EXP9_RELEASE", "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/SLEAPanalyzer_v2/releases/2026-11-05_exp9_all-batches_canonical-landmark-repair")
CORNERS <- c("tl", "tr", "br", "bl")
FPS <- 30                 # the rate the scoring assumes for one tracking row
FPS_VIDEO <- 30000 / 1001 # the rate the cameras recorded
EPM <- load_assay_config(file.path(PROJ, "config", "epm_all.yaml"), "EPM")
NOR <- load_assay_config(file.path(PROJ, "config", "nor_all.yaml"), "NOR")

# what each test should look like, and the landmarks its scores use
EXPECT <- data.table(assay = c("SocP", "NOR", "OFT", "EPM"), duration_s = c(600, 300, 600, 300),
                     width_cm = c(44, 49, 49, NA), height_cm = c(24, 49, 49, NA))
USED <- list(SocP = c("nose", "bodycentre"), NOR = c("nose", "bodycentre"), OFT = "bodycentre",
             EPM = c("bodycentre", "headcentre", "neck"))
# flags: a value outside these limits is listed for a look
LIMIT <- list(length_s = 2, last_bin = 0.9, duplicate_cm = 0.15, duplicate_lag = 1800, used_observed = 0.95,
              aspect = 0.03, arm_cm = c(4.5, 5.5), drift_cm = 1, corner_off = 0.25, camera_shift_cm = 0.5,
              moving_speed = 10, mouse = 0.2, forward = 0.8, spike_cm = 1, spike_path_pct = 1.5,
              objects_apart = 0.5, object_on_point = 0.8, budget = 0.05, oft_dropped_s = 15, contact_s = 5, z = 3.5)

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
# EXP9_LIGHT_ONLY (a regular expression on "test phase batch code") checks a subset and writes to a temporary folder
ONLY <- Sys.getenv("EXP9_LIGHT_ONLY")
if (nzchar(ONLY)) recordings <- recordings[grepl(ONLY, paste(assay, phase, batch, Code))]
OUT <- if (nzchar(ONLY)) Sys.getenv("EXP9_LIGHT_OUT", tempdir()) else RES
recordings[, md5 := unname(tools::md5sum(file))]
cat(nrow(recordings), "recordings\n")

xy <- function(tr, p) if (p %in% names(tr$data)) cbind(tr$data[[p]]$x, tr$data[[p]]$y) else matrix(NA_real_, nrow(tr$data[[1]]), 2)
observed <- function(m) is.finite(m[, 1]) & is.finite(m[, 2])
longest_run <- function(x) { r <- rle(x %in% TRUE); if (any(r$values)) max(r$lengths[r$values]) else 0L }
point_distance <- function(a, b) sqrt(sum((unlist(a) - unlist(b))^2))
norm2 <- function(m) sqrt(rowSums(m^2))
# distance of each frame from the midpoint of the frames either side
off_midpoint <- function(m) { n <- nrow(m); d <- rep(NA_real_, n); i <- 2:(n - 1)
  d[i] <- norm2(m[i, , drop = FALSE] - (m[i - 1, , drop = FALSE] + m[i + 1, , drop = FALSE]) / 2); d }
path_length <- function(m) sum(norm2(diff(m)), na.rm = TRUE)

check_one <- function(r) {
  tr <- read_tracking_csv(file = r$file, fps = FPS)
  if (r$assay == "EPM") {
    tr <- CalibrateTrackingData(tr, method = "distance", in.metric = EPM$calibration_distance_cm, points = EPM$calibration_points)
  } else {
    e <- EXPECT[assay == r$assay]
    tr <- CalibrateTrackingData(tr, method = "area", in.metric = e$width_cm * e$height_cm, points = CORNERS)
  }
  n <- nrow(tr$data$bodycentre)
  body <- xy(tr, "bodycentre"); nose <- xy(tr, "nose"); tail <- xy(tr, "tailbase")
  landmarks <- setdiff(names(tr$data), c(CORNERS, "socl", "socr", "objL", "objR", "ctl", "ctr", "cbl", "cbr", "lt", "lb", "rt", "rb"))
  share <- vapply(landmarks, function(p) mean(observed(xy(tr, p))), 0)
  used <- USED[[r$assay]]
  used_share <- vapply(used, function(p) if (p %in% names(share)) share[[p]] else 0, 0)
  out <- list(frames = n, duration_s = n / FPS_VIDEO,
              last_bin_share = if (r$assay == "OFT" && n %% 1800 > 0) (n %% 1800) / 1800 else NA_real_,
              body_observed = mean(observed(body)), nose_observed = mean(observed(nose)),
              body_longest_gap_s = longest_run(!observed(body)) / FPS,
              used_observed_min = min(used_share), used_worst = names(which.min(used_share)),
              empty_landmarks = paste(names(share)[share == 0], collapse = " "),
              cm_per_px = tr$px.to.cm)

  # arena
  cxy <- lapply(CORNERS, function(p) xy(tr, p)); names(cxy) <- CORNERS
  med <- lapply(cxy, function(m) apply(m, 2, stats::median, na.rm = TRUE))
  offs <- lapply(CORNERS, function(p) sweep(cxy[[p]], 2, med[[p]]))
  dist <- sapply(offs, norm2)
  out$corner_drift_cm <- max(apply(dist, 2, stats::quantile, 0.95, na.rm = TRUE))
  out$corner_off_share <- max(colMeans(dist > 2, na.rm = TRUE))
  bins <- (seq_len(n) - 1) %/% (10 * FPS)
  shift <- vapply(split(seq_len(n), bins), function(i) {
    dx <- stats::median(sapply(offs, function(o) stats::median(o[i, 1], na.rm = TRUE)))
    dy <- stats::median(sapply(offs, function(o) stats::median(o[i, 2], na.rm = TRUE)))
    sqrt(dx^2 + dy^2) }, 0)
  out$camera_shift_cm <- max(shift, na.rm = TRUE)
  if (r$assay == "EPM") {
    out$arm_width_cm <- point_distance(med$tl, med$tr)   # tl-bl is the 60 cm calibration span along this arm
  } else {
    out$aspect_ratio <- (point_distance(med$tl, med$tr) + point_distance(med$bl, med$br)) /
                        (point_distance(med$tl, med$bl) + point_distance(med$tr, med$br))
  }

  # the mouse while it moves
  if (any(observed(nose)) && any(observed(tail))) {
    len <- norm2(nose - tail)
    step <- rbind(NA, diff(body)); speed <- norm2(step) * FPS
    moving <- is.finite(speed) & speed > LIMIT$moving_speed
    ahead <- rowSums(step * (nose - tail))   # > 0: the body centre moves toward the nose
    out$mouse_length_cm <- stats::median(len, na.rm = TRUE)
    out$moving_length_cm <- if (sum(moving) >= 30) stats::median(len[moving], na.rm = TRUE) else NA_real_
    out$moving_forward_share <- if (sum(moving) >= 30) mean(ahead[moving] > 0, na.rm = TRUE) else NA_real_
    # a pose from another frame moves the whole skeleton off its path and back
    dn <- off_midpoint(nose); db <- off_midpoint(body); dt <- off_midpoint(tail)
    spike <- dn > LIMIT$spike_cm & db > LIMIT$spike_cm & dt > LIMIT$spike_cm
    spike[is.na(spike)] <- FALSE
    smoothed <- body
    i <- which(spike); if (length(i)) smoothed[i, ] <- (body[i - 1, , drop = FALSE] + body[i + 1, , drop = FALSE]) / 2
    path <- path_length(body)
    out$spike_rows <- sum(spike)
    out$spike_path_pct <- if (path > 0) 100 * (path - path_length(smoothed)) / path else 0
    out$body_flick_rows <- sum(db > LIMIT$spike_cm & dn < 0.5 & dt < 0.5, na.rm = TRUE)
  }

  # novel object positions
  if (r$assay == "NOR") {
    st <- tryCatch(hold_static_landmarks(tr, c("objL", "objR"), NOR$static_objects_min_separation_cm), error = function(e) NULL)
    if (!is.null(st)) {
      fixed <- list(objL = unlist(st$data$objL[1, c("x", "y")]), objR = unlist(st$data$objR[1, c("x", "y")]))
      out$object_separation_cm <- point_distance(fixed$objL, fixed$objR)
      out$objects_apart <- attr(st, "static_landmarks")$share_of_frames_used
      on_point <- vapply(names(fixed), function(p) mean(norm2(sweep(xy(tr, p), 2, fixed[[p]])) <= 2, na.rm = TRUE), 0)
      out$object_on_point_min <- min(on_point)
    }
  }
  list(values = out, body = body)
}

t0 <- Sys.time()
rows <- vector("list", nrow(recordings)); tracks <- vector("list", nrow(recordings))
for (i in seq_len(nrow(recordings))) {
  r <- recordings[i]
  v <- tryCatch(check_one(r), error = function(e) list(values = list(error = conditionMessage(e)), body = NULL))
  rows[[i]] <- cbind(r[, .(assay, phase, batch, Code, file = basename(file))], as.data.table(v$values))
  tracks[[i]] <- v$body
  if (i %% 50 == 0) cat(i, "of", nrow(recordings), format(round(Sys.time() - t0)), "\n")
}
chk <- rbindlist(rows, fill = TRUE)
if ("error" %in% names(chk)) print(chk[!is.na(error), .(assay, phase, batch, Code, error)])
for (cn in c("arm_width_cm", "aspect_ratio", "mouse_length_cm", "moving_length_cm", "moving_forward_share", "spike_rows",
             "spike_path_pct", "body_flick_rows", "object_separation_cm", "objects_apart", "object_on_point_min"))
  if (!cn %in% names(chk)) chk[, (cn) := NA_real_]   # a subset run may lack the tests that fill them

# ---- the same trial twice ------------------------------------------------------------------------------------------------
# Identical files anywhere, and within each test and batch every pair of recordings compared at the offset (within 60 s)
# where their body-centre tracks correlate best: the median distance there is about 0 cm only if they are the same trial.
dups <- recordings[, .(n = .N, files = paste(paste(assay, phase, batch, Code), collapse = " = ")), by = md5][n > 1]
dups <- if (nrow(dups)) dups[, .(pair = files, offset_frames = 0L, median_distance_cm = 0, how = "identical file")] else NULL
# cross-correlation by FFT: element k + 1 of the inverse transform is sum over t of a(t + k) b(t)
prep <- function(m) { for (k in 1:2) { v <- m[, k]; v[!is.finite(v)] <- stats::median(v, na.rm = TRUE); m[, k] <- v - mean(v) }; m }
median_gap <- function(a, b, lag) {   # b(t) against a(t + lag)
  ia <- seq_len(nrow(a)); ib <- ia - lag; ok <- ib >= 1 & ib <= nrow(b)
  if (sum(ok) < 300) return(Inf)
  stats::median(norm2(a[ia[ok], , drop = FALSE] - b[ib[ok], , drop = FALSE]), na.rm = TRUE)
}
lag_pairs <- list()
for (g in split(seq_len(nrow(recordings)), recordings[, paste(assay, batch)])) {
  g <- g[!vapply(tracks[g], is.null, TRUE)]
  if (length(g) < 2) next
  # one transform per recording; x and y correlations summed per pair
  N <- stats::nextn(2 * max(vapply(tracks[g], nrow, 0L)))
  Fx <- lapply(tracks[g], function(m) { p <- prep(m); list(x = stats::fft(c(p[, 1], rep(0, N - nrow(p)))), y = stats::fft(c(p[, 2], rep(0, N - nrow(p))))) })
  lags <- -LIMIT$duplicate_lag:LIMIT$duplicate_lag; at <- ifelse(lags >= 0, lags + 1, N + lags + 1)
  for (a in seq_along(g)[-length(g)]) for (b in (a + 1):length(g)) {
    cc <- Re(stats::fft(Fx[[a]]$x * Conj(Fx[[b]]$x) + Fx[[a]]$y * Conj(Fx[[b]]$y), inverse = TRUE))[at]
    # the correlation peak of smooth tracks can sit a frame or two off the true offset: take the closest fit around it
    near <- lags[which.max(cc)] + (-3:3)
    gaps <- vapply(near, function(l) median_gap(tracks[[g[a]]], tracks[[g[b]]], l), 0)
    lag <- near[which.min(gaps)]; d <- min(gaps)
    if (d < LIMIT$duplicate_cm) lag_pairs[[length(lag_pairs) + 1]] <- data.table(
      pair = paste(recordings[g[a], paste(assay, phase, batch, Code)], "=", recordings[g[b], paste(assay, phase, batch, Code)]),
      offset_frames = lag, median_distance_cm = d, how = "same track at an offset")
  }
}
duplicates <- rbindlist(c(list(dups), lag_pairs), fill = TRUE)
if (!nrow(duplicates)) duplicates <- data.table(pair = character(), offset_frames = integer(), median_distance_cm = numeric(), how = character())
fwrite(duplicates, file.path(OUT, "light_checks_duplicates.csv"))
dup_members <- unique(unlist(strsplit(duplicates$pair, " = ")))
chk[, duplicate := paste(assay, phase, batch, Code) %in% dup_members]

# ---- the release: animal metadata, readouts, plus maze zone budget, open field dropped time ----------------------------
animals <- fread(file.path(RELEASE, "data", "exp9_sleap_all_animals.tsv"))
chk <- merge(chk, animals[, .(Code, Sex)], by = "Code", all.x = TRUE)
budget <- rbindlist(lapply(list.files(file.path(RELEASE, "assay_summaries", "EPM"), "^Report[.]csv$", recursive = TRUE, full.names = TRUE),
  function(f) { x <- fread(f); x[, .(Code = str_extract(file, "^[A-Za-z0-9]{4}"),
    no_zone_share = 1 - (bodycentre.open.total.time + bodycentre.closed.total.time + bodycentre.center.total.time) / bodycentre.valid.time)] }))
chk <- merge(chk, budget[, .(Code, assay = "EPM", no_zone_share)], by = c("Code", "assay"), all.x = TRUE)
dropped <- rbindlist(lapply(list.files(file.path(RELEASE, "assay_summaries", "OFT"), "^combined_qc[.]csv$", recursive = TRUE, full.names = TRUE),
  function(f) fread(f)[stage == "cleaned" & point == "bodycentre", .(Code = file, oft_dropped_s = missing_fraction * duration_s)]))
chk <- merge(chk, dropped[, .(Code, assay = "OFT", oft_dropped_s)], by = c("Code", "assay"), all.x = TRUE)
chk <- merge(chk, animals[, .(Code, assay = "NOR", phase = "NOV", total_contact_s = sleap_NOR_total_s)], by = c("Code", "assay", "phase"), all.x = TRUE)
chk <- merge(chk, animals[, .(Code, assay = "SocP", phase = "S1", s1_pref = sleap_SocP_S1_pref_index)], by = c("Code", "assay", "phase"), all.x = TRUE)

# ---- frame-level checks from the plausibility screen (information) -----------------------------------------------------
scr <- fread(file.path(RES, "plausibility_summary.csv"))
chk <- merge(chk, scr[, .(assay, phase, batch, Code, jump_s, skeleton_s, outside_arena_s, raised_s)],
             by = c("assay", "phase", "batch", "Code"), all.x = TRUE)

# ---- readout outliers: robust z within test and sex, on the trial they come from --------------------------------------
readouts <- data.table(assay = c("EPM", "EPM", "OFT", "OFT", "NOR", "NOR", "SocP", "SocP"),
                       phase = c("EPM", "EPM", "OFT", "OFT", "NOV", "NOV", "S1", "S2"),
                       column = c("sleap_EPM_open_frac", "sleap_EPM_distance_cm", "sleap_OFT_center_pct", "sleap_OFT_distance_cm",
                                  "sleap_NOR_D2", "sleap_NOR_distance_cm", "sleap_SocP_S1_pref_index", "sleap_SocP_S2_pref_index"))
rz <- rbindlist(lapply(seq_len(nrow(readouts)), function(k) {
  x <- animals[, .(Code, Sex, value = get(readouts$column[k]))]
  x[, z := (value - stats::median(value, na.rm = TRUE)) / stats::mad(value, na.rm = TRUE), by = Sex]
  x[is.finite(z) & abs(z) > LIMIT$z, .(Code, assay = readouts$assay[k], phase = readouts$phase[k],
                                       readout = sub("^sleap_", "", readouts$column[k]), z = round(z, 1))] }))
rz <- rz[, .(readout_outliers = paste0(readout, " z ", z, collapse = "; ")), by = .(Code, assay, phase)]
chk <- merge(chk, rz, by = c("Code", "assay", "phase"), all.x = TRUE)

# ---- the video QC, when it has been run ----------------------------------------------------------------------------------
# tools/video_qc run on the manifest of 32_video_qc_manifest.R with out_dir results/video_qc (the command is in 32). Its
# flags join as video:<check>:<flag>; video_flags is "none" for a recording checked without a flag, empty for one not
# checked, and video_checks names the checks that ran on it.
VIDEO_QC <- Sys.getenv("EXP9_VIDEO_QC", file.path(RES, "video_qc", "video_qc.csv"))
vq <- if (file.exists(VIDEO_QC) && file.size(VIDEO_QC) > 2) tryCatch(fread(VIDEO_QC, colClasses = "character"), error = function(e) NULL)
if (file.exists(VIDEO_QC) && !all(c("id", "flags") %in% names(vq))) {
  cat("video QC file empty or not readable, video flags not joined:", VIDEO_QC, "\n"); vq <- NULL
} else if (!file.exists(VIDEO_QC) && length(list.files(dirname(VIDEO_QC), "^video_qc_.*[.]csv$")))
  cat("video QC tables but no video_qc.csv in", dirname(VIDEO_QC), "(a run in progress or stopped): video flags not joined\n")
if (!is.null(vq)) {
  if (!"checks" %in% names(vq)) vq[, checks := NA_character_]
  vq[flags == "", flags := "none"]
  m <- match(paste(chk$assay, chk$phase, chk$batch, chk$Code), vq$id)
  chk[, `:=`(video_checks = vq$checks[m], video_flags = vq$flags[m])]
  cat(sum(!is.na(chk$video_flags)), "of", nrow(chk), "recordings have video QC results in", VIDEO_QC, "\n")
  # the oldest table it combines decides whether the video QC is older than the tracking
  tabs <- file.path(dirname(VIDEO_QC), paste0("video_qc_", unique(unlist(strsplit(stats::na.omit(vq$checks), ";", fixed = TRUE))), ".csv"))
  if (any(file.mtime(recordings$file) > min(file.mtime(c(VIDEO_QC, tabs[file.exists(tabs)])))))
    cat("  tracking files newer than the video QC: rerun 32_video_qc_manifest.R and tools/video_qc\n")
}

# ---- flags ---------------------------------------------------------------------------------------------------------------
chk <- merge(chk, EXPECT[, .(assay, expected_s = duration_s, expected_aspect = width_cm / height_cm)], by = "assay")
chk[, moving_length_rel := moving_length_cm / stats::median(moving_length_cm, na.rm = TRUE), by = .(assay, Sex, batch)]
flag <- function(cond, name) ifelse(cond %in% TRUE, name, NA_character_)
chk[, flags := apply(cbind(
  flag(abs(duration_s - expected_s) > LIMIT$length_s, "length"),
  flag(last_bin_share < LIMIT$last_bin, "partial_last_bin"),
  flag(duplicate, "duplicate"),
  flag(used_observed_min < LIMIT$used_observed, "landmark_missing"),
  flag(abs(aspect_ratio / expected_aspect - 1) > LIMIT$aspect, "arena_shape"),
  flag(arm_width_cm < LIMIT$arm_cm[1] | arm_width_cm > LIMIT$arm_cm[2], "arm_width"),
  flag(corner_drift_cm > LIMIT$drift_cm | corner_off_share > LIMIT$corner_off, "corner"),
  flag(camera_shift_cm > LIMIT$camera_shift_cm, "camera_shift"),
  flag(abs(moving_length_rel - 1) > LIMIT$mouse, "mouse_size"),
  flag(moving_forward_share < LIMIT$forward, "nose_tail_order"),
  flag(spike_path_pct > LIMIT$spike_path_pct, "pose_spikes"),
  flag(objects_apart < LIMIT$objects_apart | object_on_point_min < LIMIT$object_on_point, "objects"),
  flag(no_zone_share > LIMIT$budget, "off_maze"),
  flag(oft_dropped_s > LIMIT$oft_dropped_s, "oft_dropped"),
  flag(total_contact_s < LIMIT$contact_s, "low_contact"),
  flag(s1_pref < 0, "s1_preference_below_0"),
  flag(!is.na(readout_outliers), "readout_outlier")), 1, function(x) paste(stats::na.omit(x), collapse = "; "))]
if ("video_flags" %in% names(chk))
  chk[, flags := mapply(function(f, v) paste(c(f[nzchar(f)], if (!v %in% c(NA, "none")) paste0("video:", strsplit(v, ";", fixed = TRUE)[[1]])),
                                             collapse = "; "), flags, video_flags, USE.NAMES = FALSE)]
setcolorder(chk, c("assay", "phase", "batch", "Code", "Sex", "flags"))
setorder(chk, assay, phase, batch, Code)
num <- names(chk)[vapply(chk, is.double, TRUE)]
chk[, (num) := lapply(.SD, signif, 4), .SDcols = num]
fwrite(chk, file.path(OUT, "light_checks_recordings.csv"))

own <- function(x) gsub("video:[^;]*(; )?", "", x)   # this script's own flags; the video ones are counted apart
flagged <- function(x, name) sum(grepl(name, own(x), fixed = TRUE))
rng <- function(x, d = 2) if (all(is.na(x))) NA_character_ else paste(round(range(x, na.rm = TRUE), d), collapse = "-")
summ <- chk[, .(recordings = .N, flagged = sum(nzchar(trimws(sub("; $", "", own(flags))))),
  length = flagged(flags, "length"), duration_range_s = rng(duration_s, 1), partial_last_bin = flagged(flags, "partial_last_bin"),
  duplicate = flagged(flags, "duplicate"), landmark_missing = flagged(flags, "landmark_missing"),
  used_observed_min_pct = round(100 * min(used_observed_min, na.rm = TRUE), 2),
  aspect_dev_max_pct = if (all(is.na(aspect_ratio))) NA_real_ else round(100 * max(abs(aspect_ratio / expected_aspect - 1), na.rm = TRUE), 2),
  arm_width_range_cm = rng(arm_width_cm), corner_drift_max_cm = round(max(corner_drift_cm, na.rm = TRUE), 2),
  corner_off_share_max = round(max(corner_off_share, na.rm = TRUE), 3), camera_shift_max_cm = round(max(camera_shift_cm, na.rm = TRUE), 2),
  corner = flagged(flags, "corner"), camera_shift = flagged(flags, "camera_shift"),
  moving_length_median_cm = round(stats::median(moving_length_cm, na.rm = TRUE), 2), mouse_size = flagged(flags, "mouse_size"),
  forward_share_min = round(min(moving_forward_share, na.rm = TRUE), 3), nose_tail_order = flagged(flags, "nose_tail_order"),
  spike_path_max_pct = round(max(spike_path_pct, na.rm = TRUE), 2), pose_spikes = flagged(flags, "pose_spikes"),
  objects_apart_min = if (all(is.na(objects_apart))) NA_real_ else round(min(objects_apart, na.rm = TRUE), 3),
  object_on_point_min = if (all(is.na(object_on_point_min))) NA_real_ else round(min(object_on_point_min, na.rm = TRUE), 3),
  objects = flagged(flags, "objects"),
  off_maze_max_pct = if (all(is.na(no_zone_share))) NA_real_ else round(100 * max(no_zone_share, na.rm = TRUE), 2),
  oft_dropped_max_s = if (all(is.na(oft_dropped_s))) NA_real_ else round(max(oft_dropped_s, na.rm = TRUE), 1),
  oft_dropped = flagged(flags, "oft_dropped"), low_contact = flagged(flags, "low_contact"),
  readout_outliers = flagged(flags, "readout_outlier")), by = .(assay, phase)]
if ("video_flags" %in% names(chk))
  summ <- merge(summ, chk[, .(video_checked = sum(!is.na(video_flags)), video = sum(grepl("video:", flags, fixed = TRUE))), by = .(assay, phase)],
                by = c("assay", "phase"), sort = FALSE)
fwrite(summ, file.path(OUT, "light_checks_summary.csv"))
print(as.data.frame(summ))
cat("\nsame trial twice:\n"); print(as.data.frame(duplicates))
cat("\nflagged recordings:\n"); print(as.data.frame(chk[flags != "", .(assay, phase, batch, Code, flags)]))
