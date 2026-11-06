# ============================================================================
# 26_frame_agreement.R
#
# The automated social-preference contact against the manual scoring of batch
# 1, frame by frame. The raw BORIS project holds the events (interaction,
# rearing, climbing per cup) of 18 of the 20 animals: camera 2 (10 animals,
# where the head length restoration was developed) and camera 1 (8 animals,
# held out), each in the habituation trial and in phases 1 and 2. Manual
# contact = interaction + rearing at a cup; climbing is not contact.
#
#   Rscript 26_frame_agreement.R <reference_release>
#
# reference_release: the SLEAPanalyzer release the automated values without
# the restoration must reproduce (its assay_summaries/SocP/B1).
#
# Each recording is processed as DLCA_SocP does (calibration, landmark repair
# with its head check, interpolation, on-cup detection; without the head check it
# must reproduce the reference release) and contact is scored with three noses:
# without restoration, restored as configured (restore_head_length(), heads
# below 0.75 of their median length), and restored in every head shorter than
# its median (sensitivity). BORIS times are seconds of video at 30 frames per
# second, like the tracked trial. Each observation is aligned to its
# recording by the offset (within 150 s) that best correlates its manual
# contact with the ear midpoint being within 7.5 cm of the cup, a signal that
# moves with the head but depends on neither nose (the restoration moves only
# the nose); only the frames the observation covers are compared. Two phase-1 observations were scored on a later chapter of the
# video and cover a window shifted against the trial (C7E2, A1L9). OQ752
# (A1L9) phase 2 is left out: its record describes its phase-1 trial
# (manual_exclusions.csv of the validation folder). One habituation trial
# (C7E2) was scored twice, which gives one manual-against-manual comparison.
#
# Writes results/socp_frame_*.csv and results/socp_bout_*.csv:
#   socp_frame_offsets.csv        per observation: offset, correlation at it, overlap
#   socp_frame_summary.csv        frame agreement (precision, recall, F1 with 95% intervals over animals) and totals (r, CCC, bias,
#                                 limits of agreement), by nose, camera and cup
#   socp_frame_totals.csv         per cup session: manual and automated contact over the compared frames, each nose
#   socp_frame_sweep.csv          frame agreement and totals against the contact distance, by nose, camera and cup
#   socp_frame_by_side.csv        manual contact by side of the cup: seconds, nose distance, scored by both / missed / added
#   socp_frame_lid_edge.csv       the place along the cup lid: visible head length, scored and unscored
#   socp_frame_distance.csv       nose distance histogram, frames with and without manual contact
#   socp_frame_positions.csv      nose positions around the cup (0.5 cm grid; x toward the other cup, y toward the lid side)
#   socp_frame_human_pair.csv     the habituation trial scored twice: manual against manual and each against the automated
#   socp_bout_sessions.csv        per cup session: manual and automated bouts, mean bout and latency
#   socp_bout_agreement.csv       their agreement against the automated merge gap and minimum bout
# ============================================================================

suppressMessages({ library(sp); library(stringr); library(yaml); library(data.table); library(jsonlite) })
SCRIPTS <- local({ a <- commandArgs(trailingOnly = FALSE); p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/") })
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
REPO <- normalizePath(file.path(PROJ, "..", ".."), winslash = "/")
Sys.setenv(SLEAP_ANALYZER_REPO_ROOT = REPO)
source(file.path(REPO, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "Behavioral_Metrics_Phase1.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "core", "assay_config.R"))
args <- commandArgs(TRUE); if (length(args) != 1) stop("usage: Rscript 26_frame_agreement.R <reference_release>")
REL <- args[1]
STAGE <- file.path(Sys.getenv("EXP9_STAGING", "C:/Users/topohl/Documents/exp9_sleap_staging"), "SocP", "B1")
E9 <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress"
BORIS <- file.path(E9, "Raw Data/Behavior/B1/SocP/BORIS/SocPref_new.boris")
NOVEL <- file.path(E9, "Raw Data/Behavior/B1/SocP/novelLoc%s.txt")
OUT <- file.path(PROJ, "results")
CFG <- list(S1 = "socp_all.yaml", S2 = "socp_all.yaml", HAB = "socp_hab_all.yaml")
configs <- lapply(CFG, function(f) load_assay_config(file.path(PROJ, "config", f), "SocP"))
stopifnot("the configurations restore the head length" = all(vapply(configs, function(x) isTRUE(x$head_length_restore), TRUE)))
config <- configs$S1
FPS <- config$fps; stopifnot(FPS == 30)
CAMERA2 <- c("A1L9", "C7E2", "F9L3", "K4T2", "M1M3", "P9X0", "R8O3", "T2H7", "V2P0", "Y0F3")
MAX_LAG_S <- 150; EAR_NEAR_CM <- 7.5; MIN_PEAK_R <- 0.25
NOSES <- c("without restoration", "restored", "restored, every short head")
dist <- function(x1, y1, x2, y2) sqrt((x1 - x2)^2 + (y1 - y2)^2)
ccc <- function(x, y) { n <- length(x); sxy <- stats::cov(x, y) * (n - 1) / n; sx <- stats::var(x) * (n - 1) / n; sy <- stats::var(y) * (n - 1) / n
  2 * sxy / (sx + sy + (mean(x) - mean(y))^2) }   # Lin's population form, as 05, 12, 24

# ---------------------------------------------------------------- manual observations
obs <- fromJSON(BORIS, simplifyVector = FALSE)$observations
O <- rbindlist(lapply(names(obs), function(k) {
  e <- obs[[k]]$events; if (!length(e)) return(NULL)
  subj <- unique(vapply(e, function(z) z[[2]], "")); stopifnot(length(subj) == 1)
  phase <- if (grepl("_hab$", k)) "HAB" else if (grepl("_stim1$", k)) "S1" else if (grepl("_stim2$", k)) "S2" else if (k == "GX0787_1.mp4") "HAB" else NA
  data.table(observation = k, code = subj, phase = phase, second_scoring = k == "GX0787_1.mp4") }))
O <- O[!is.na(phase) & !(code == "A1L9" & phase == "S2")]
O[, camera := fifelse(code %in% CAMERA2, "camera 2", "camera 1")]
manual_events <- function(k) {
  ob <- obs[[k]]; iv <- unlist(ob[["observation time interval"]]); e <- ob$events
  t <- vapply(e, function(z) as.numeric(z[[1]]), 0) - iv[1]; b <- vapply(e, function(z) z[[3]], "")
  n <- ceiling(round((iv[2] - iv[1]) * FPS, 6)); fr <- seq_len(n) - 1
  state <- function(beh) { k <- which(b == beh); stopifnot("state events come in start/stop pairs" = length(k) %% 2 == 0)
    s <- t[k[c(TRUE, FALSE)]]; f <- t[k[c(FALSE, TRUE)]]; x <- logical(n)
    for (i in seq_along(s)) x[fr >= round(s[i] * FPS) & fr < round(f[i] * FPS)] <- TRUE; x }
  list(L = state("l interact") | state("l rear"), R = state("r interact") | state("r rear"), rearL = state("l rear"), rearR = state("r rear"), n = n)
}

# ---------------------------------------------------------------- per frame, as DLCA_SocP
qcDecisions <- read_qc_decisions(config$qc_decisions_file)
process <- function(code, ph) {
  cf <- configs[[ph]]; name <- paste0(code, "_", ph)
  tr <- read_tracking_csv(file = file.path(STAGE, ph, paste0(name, ".csv")), fps = cf$fps)
  raw <- CalibrateTrackingData(tr, method = "area", in.metric = cf$arena_width_cm * cf$arena_height_cm, points = cf$arena_corner_names)
  prepare <- function(head_check) {
    t <- repair_landmarks(raw, landmarks = cf$landmark_repair_landmarks, targets = cf$landmark_repair_targets,
      margin = cf$landmark_repair_margin_cm, margin_fraction = cf$landmark_repair_margin_fraction, min_share = cf$landmark_repair_min_share,
      min_landmarks = cf$landmark_repair_min_landmarks, max_anchor_error = cf$landmark_repair_max_anchor_error_cm,
      corners = cf$arena_corner_names, max_outside = cf$landmark_repair_outside_cm, head_check = head_check)
    t <- interpolate_tracking(t, landmarks = c("nose", "bodycentre"), max_gap_s = cf$max_interpolation_gap_s, likelihood_cutoff = cf$likelihood_cutoff)
    size <- apparent_body_size(t)
    on <- lapply(c(left = "socl", right = "socr"), function(cup) {
      f <- detect_on_object(t, cup, cf$arena_corner_names, footprint_radius = cf$on_cup_footprint_radius_cm, beyond = cf$on_cup_beyond_cm,
        max_distance = cf$on_cup_max_distance_cm, min_apparent_size = cf$on_cup_min_apparent_size, merge_gap_s = cf$on_cup_merge_gap_s,
        min_episode_s = cf$on_cup_min_episode_s, size_window_s = cf$on_cup_size_window_s, size = size)
      apply_qc_decisions(f$on_object, qcDecisions, name, cup, cf$fps) })
    list(t = t, on = on)
  }
  score <- function(t, on) compute_socp_metrics(t, "R", cf$fps, contact_distance = cf$contact_distance_cm,
    body_exclusion_distance = cf$body_exclusion_distance_cm, proximity_range = cf$proximity_range_cm, threshold_unit = "cm", on_cup = on)$summary
  released <- prepare(FALSE); x <- prepare(TRUE); tr <- x$t; on <- x$on
  variants <- list(tr, restore_head_length(tr, below = cf$head_length_restore_below), restore_head_length(tr, below = 1))
  names(variants) <- NOSES
  sums <- lapply(variants, score, on = on)
  g <- function(t, p, xy) t$data[[p]][[xy]]
  d <- data.table(code = code, phase = ph, frame = seq_along(g(tr, "nose", "x")), valid = landmark_validity(tr, c("nose", "bodycentre", "socl", "socr")),
                  on_l = on$left, on_r = on$right, body_x = g(tr, "bodycentre", "x"), body_y = g(tr, "bodycentre", "y"),
                  socl_x = g(tr, "socl", "x"), socl_y = g(tr, "socl", "y"), socr_x = g(tr, "socr", "x"), socr_y = g(tr, "socr", "y"),
                  em_x = (g(tr, "leftEar", "x") + g(tr, "rightEar", "x")) / 2, em_y = (g(tr, "leftEar", "y") + g(tr, "rightEar", "y")) / 2,
                  head_ratio = variants[[2]]$head_length$length / variants[[2]]$head_length$median)
  for (v in seq_along(NOSES)) { d[, (paste0("nx", v)) := g(variants[[v]], "nose", "x")]; d[, (paste0("ny", v)) := g(variants[[v]], "nose", "y")] }
  rs <- score(released$t, released$on)
  attr(d, "totals") <- data.table(nose = "as released", left = rs$contactLeft, right = rs$contactRight)
  d
}
REC <- unique(O[, .(code, phase, camera)])
frames <- list(); totals_check <- list()
for (i in seq_len(nrow(REC))) {
  d <- process(REC$code[i], REC$phase[i]); frames[[i]] <- d
  totals_check[[i]] <- cbind(REC[i], attr(d, "totals"))
}
totals_check <- rbindlist(totals_check)
rel <- rbindlist(lapply(c("HAB", "S1", "S2"), function(ph) fread(file.path(REL, "assay_summaries", "SocP", "B1", ph, "combined_output.csv"))))
chk <- totals_check[, file := paste0(code, "_", phase)][rel, on = "file", nomatch = 0]
stopifnot("without the head check and the restoration the release values are reproduced" = nrow(chk) == nrow(REC),
          max(abs(chk$left - chk$contactLeft), abs(chk$right - chk$contactRight)) < 1e-6)
cat(nrow(REC), "recordings processed; the release values are reproduced without the head check and the restoration\n")

# ---------------------------------------------------------------- alignment (rule-independent signal)
xcor <- function(x, y, max_lag) {   # correlation of x (manual frames) with y shifted: tracking row = manual frame + lag
  x <- x - mean(x); y <- y - mean(y); nx <- length(x); ny <- length(y); N <- 2^ceiling(log2(nx + ny + 2 * max_lag))
  cc <- Re(stats::fft(Conj(stats::fft(c(x, rep(0, N - nx)))) * stats::fft(c(y, rep(0, N - ny))), inverse = TRUE)) / N
  lags <- -max_lag:max_lag; v <- cc[ifelse(lags >= 0, lags + 1, N + lags + 1)]
  ov <- pmin(nx, ny - lags) - pmax(0, -lags)
  list(lags = lags, r = v / (pmax(ov, 1) * stats::sd(x) * stats::sd(y)))
}
nearsig <- function(d, side) { s <- if (side == "L") "socl" else "socr"
  z <- dist(d$em_x, d$em_y, d[[paste0(s, "_x")]], d[[paste0(s, "_y")]]) <= EAR_NEAR_CM; z[is.na(z)] <- FALSE; as.numeric(z) }
M <- list(); offsets <- list(); human <- list()
for (j in seq_len(nrow(O))) {
  o <- O[j]; d <- frames[[which(REC$code == o$code & REC$phase == o$phase)]]; me <- manual_events(o$observation)
  a <- xcor(as.numeric(me$L), nearsig(d, "L"), MAX_LAG_S * FPS); b <- xcor(as.numeric(me$R), nearsig(d, "R"), MAX_LAG_S * FPS)
  r <- (a$r + b$r) / 2; k <- which.max(r); lag <- a$lags[k]
  away <- abs(a$lags - lag) > 10 * FPS                     # a rival peak is a different episode, not the slope of this one
  offsets[[j]] <- data.table(o[, .(observation, code, phase, camera, second_scoring)], offset_s = lag / FPS, r_at_offset = r[k],
                             r_next_best = max(r[away]), r_at_zero = r[a$lags == 0], compared_s = NA_real_)
  rows <- seq_len(nrow(d)); m <- rows - 1 - lag + 1; inside <- m >= 1 & m <= me$n
  offsets[[j]]$compared_s <- sum(inside) / FPS
  long <- rbind(d[inside, .(code, phase, frame, cup = "L", m = me$L[m[inside]], rear = me$rearL[m[inside]], valid, on = on_l,
                            sx = socl_x, sy = socl_y, ox = socr_x, oy = socr_y, body_x, body_y, head_ratio,
                            nx1, ny1, nx2, ny2, nx3, ny3)],
                d[inside, .(code, phase, frame, cup = "R", m = me$R[m[inside]], rear = me$rearR[m[inside]], valid, on = on_r,
                            sx = socr_x, sy = socr_y, ox = socl_x, oy = socl_y, body_x, body_y, head_ratio,
                            nx1, ny1, nx2, ny2, nx3, ny3)])
  long[, `:=`(observation = o$observation, camera = o$camera, full_window = abs(lag) <= 2 * FPS)]
  if (o$second_scoring) human[[length(human) + 1]] <- long else M[[length(M) + 1]] <- long
}
offsets <- rbindlist(offsets)
fwrite(offsets, file.path(OUT, "socp_frame_offsets.csv"))
print(offsets[, .(observation, code, phase, offset_s, r = round(r_at_offset, 3), next_best = round(r_next_best, 3), at_zero = round(r_at_zero, 3), compared_s)])
stopifnot("every alignment has a clear peak" = all(offsets$r_at_offset >= MIN_PEAK_R), all(offsets$r_at_offset > 1.5 * offsets$r_next_best))
cat(sprintf("%d observations aligned: offsets %s s; shifted windows: %s\n", nrow(offsets), paste(range(round(offsets$offset_s, 2)), collapse = " to "),
            paste(offsets[abs(offset_s) > 2, sprintf("%s %s %+.1f s", code, phase, offset_s)], collapse = ", ")))
M <- rbindlist(M)

# cup role and geometry: x toward the other cup, y flipped for the right cup so both lids lie at y > 0
nl <- rbindlist(lapply(c("HAB", "S1", "S2"), function(ph) fread(sprintf(NOVEL, ph))[, .(code = Code, phase = ph, NovelLoc)]))
prep <- function(D) {
  D <- nl[D, on = .(code, phase)]
  D[, role := fifelse(phase == "HAB", "habituation", fifelse((NovelLoc == "R") == (cup == "L"), "stranger cup", "other cup"))]
  D[, flip := fifelse(cup == "L", 1, -1)]
  D[, b := dist(body_x, body_y, sx, sy)]
  for (v in 1:3) { D[, (paste0("d", v)) := dist(get(paste0("nx", v)), get(paste0("ny", v)), sx, sy)] }
  D[, `:=`(rx = (nx1 - sx) * sign(ox - sx), ry = (ny1 - sy) * flip, hx = (nx2 - sx) * sign(ox - sx), hy = (ny2 - sy) * flip)]
  D[, ok := valid & !on & b > config$body_exclusion_distance_cm]; D[is.na(ok), ok := FALSE]
  for (v in 1:3) { D[, (paste0("a", v)) := ok & !is.na(get(paste0("d", v))) & get(paste0("d", v)) <= config$contact_distance_cm] }
  D
}
M <- prep(M); H <- prep(rbindlist(human))

# ---------------------------------------------------------------- frame agreement and totals
set.seed(26)
boot_f1 <- function(r, reps = 2000) { an <- unique(r$code); if (length(an) < 3) return(c(NA_real_, NA_real_))
  s <- replicate(reps, { k <- sample(an, replace = TRUE); x <- r[data.table(code = k), on = "code", allow.cartesian = TRUE]
    2 * sum(x$tp) / (2 * sum(x$tp) + sum(x$fn) + sum(x$fp)) }); stats::quantile(s, c(.025, .975)) }
per_session <- function(D, pred) D[, .(tp = sum(m & get(pred)), fn = sum(m & !get(pred)), fp = sum(!m & get(pred)),
                                       manual_s = sum(m) / FPS, auto_s = sum(get(pred)) / FPS), by = .(code, phase, cup, role, camera, full_window)]
agree <- function(r) { ci <- boot_f1(r)
  data.table(cup_sessions = nrow(r), animals = uniqueN(r$code), precision = sum(r$tp) / (sum(r$tp) + sum(r$fp)), recall = sum(r$tp) / (sum(r$tp) + sum(r$fn)),
             f1 = 2 * sum(r$tp) / (2 * sum(r$tp) + sum(r$fn) + sum(r$fp)), f1_low = ci[1], f1_high = ci[2], manual_s = sum(r$manual_s),
             missed_s = sum(r$fn) / FPS, added_s = sum(r$fp) / FPS, r = stats::cor(r$auto_s, r$manual_s), ccc = ccc(r$auto_s, r$manual_s),
             bias_s = mean(r$auto_s - r$manual_s), loa_s = 1.96 * stats::sd(r$auto_s - r$manual_s), mad_s = mean(abs(r$auto_s - r$manual_s))) }
groupings <- list(list(), list("camera"), list("camera", "role"))
summ <- rbindlist(lapply(1:3, function(v) { ps <- per_session(M, paste0("a", v))
  rbindlist(lapply(groupings, function(g) { if (!length(g)) return(cbind(nose = NOSES[v], camera = "both", role = "all", agree(ps)))
    keys <- unique(ps[, unlist(g), with = FALSE])
    rbindlist(lapply(seq_len(nrow(keys)), function(i) { r <- ps[keys[i], on = unlist(g)]
      cbind(nose = NOSES[v], camera = keys$camera[i], role = if ("role" %in% g) keys$role[i] else "all", agree(r)) })) })) }))
fwrite(summ, file.path(OUT, "socp_frame_summary.csv"))
print(summ[role == "all", .(nose, camera, cup_sessions, animals, f1 = round(f1, 3), ci = sprintf("%.3f-%.3f", f1_low, f1_high), ccc = round(ccc, 3),
                            bias_s = round(bias_s, 1), loa_s = round(loa_s, 1))])
tot <- rbindlist(lapply(1:3, function(v) cbind(nose = NOSES[v], per_session(M, paste0("a", v)))))
fwrite(tot, file.path(OUT, "socp_frame_totals.csv"))

# ---------------------------------------------------------------- contact distance sweep
sweep <- rbindlist(lapply(seq(4, 9, by = .25), function(cd) rbindlist(lapply(1:2, function(v) {
  M[, z := ok & !is.na(get(paste0("d", v))) & get(paste0("d", v)) <= cd]; ps <- per_session(M, "z")
  rbindlist(lapply(list(c("camera"), c("camera", "role")), function(g) {
    ps[, { r <- .SD; data.table(cup_sessions = .N, f1 = 2 * sum(r$tp) / (2 * sum(r$tp) + sum(r$fn) + sum(r$fp)), ccc = ccc(r$auto_s, r$manual_s),
                                bias_s = mean(r$auto_s - r$manual_s), loa_s = 1.96 * stats::sd(r$auto_s - r$manual_s)) }, by = g][, `:=`(nose = NOSES[v], contact_distance_cm = cd)] }), fill = TRUE) }))))
sweep[is.na(role), role := "all"]; M[, z := NULL]
fwrite(sweep, file.path(OUT, "socp_frame_sweep.csv"))

# ---------------------------------------------------------------- where: side of the cup, the lid edge, distances, positions
M[, angle := abs(((atan2(ry, rx) * 180 / pi) + 180) %% 360 - 180)]
M[, side := cut(angle, c(0, 45, 90, 135, 180), include.lowest = TRUE, labels = c("0-45", "45-90", "90-135", "135-180"))]
side <- M[valid == TRUE & !is.na(side), .(manual_s = sum(m) / FPS, both_s = sum(m & a1) / FPS, missed_s = sum(m & !a1) / FPS, added_s = sum(!m & a1) / FPS,
          both_restored_s = sum(m & a2) / FPS, missed_restored_s = sum(m & !a2) / FPS, added_restored_s = sum(!m & a2) / FPS,
          median_nose_cm = stats::median(d1[m], na.rm = TRUE)), by = .(camera, role, side)][order(camera, role, side)]
fwrite(side, file.path(OUT, "socp_frame_by_side.csv"))
M[, lid_edge := valid & d1 > 6.3 & d1 < 11 & ry > 3.8 & rx < -2.5]
lid <- M[lid_edge == TRUE & is.finite(head_ratio), .(time_s = .N / FPS, median_head_ratio = stats::median(head_ratio), share_head_below_75 = mean(head_ratio < .75),
                                                    contact_restored = mean(a2)), by = .(camera, manual_contact = m)][order(camera, manual_contact)]
fwrite(lid, file.path(OUT, "socp_frame_lid_edge.csv")); print(lid, digits = 3)
br <- seq(0, 20, by = .25)
hist <- rbindlist(lapply(1:2, function(v) M[valid == TRUE & get(paste0("d", v)) < 20, .(time_s = .N / FPS),
  by = .(camera, manual_contact = m, bin_cm = br[findInterval(get(paste0("d", v)), br)] + .125)][, nose := NOSES[v]]))
fwrite(hist[order(nose, camera, manual_contact, bin_cm)], file.path(OUT, "socp_frame_distance.csv"))
st <- function(m, a) fifelse(m & a, "both", fifelse(m, "manual only", fifelse(a, "automated only", "neither")))
pos <- rbind(M[valid == TRUE & d1 < 14, .(nose = NOSES[1], camera, x = round(rx * 2) / 2, y = round(ry * 2) / 2, status = st(m, a1))],
             M[valid == TRUE & d2 < 14, .(nose = NOSES[2], camera, x = round(hx * 2) / 2, y = round(hy * 2) / 2, status = st(m, a2))])
fwrite(pos[status != "neither", .(time_s = .N / FPS), by = .(nose, camera, status, x, y)], file.path(OUT, "socp_frame_positions.csv"))
cat(sprintf("manual contact beyond 14 cm of the cup point (not in the positions table): %.1f s\n", M[valid == TRUE & m & d1 >= 14, .N] / FPS))

# ---------------------------------------------------------------- the habituation trial scored twice
H1 <- M[observation == "GX010787_hab"]; H2 <- H[observation == "GX0787_1.mp4"]
hp <- merge(H1[, .(frame, cup, m1 = m, a1, a2)], H2[, .(frame, cup, m2 = m)], by = c("frame", "cup"))
f1 <- function(x, y) 2 * sum(x & y) / (2 * sum(x & y) + sum(x & !y) + sum(!x & y))
human <- data.table(comparison = c("manual 1 against manual 2", "automated (without restoration) against manual 1", "automated (without restoration) against manual 2",
                                   "automated (restored) against manual 1", "automated (restored) against manual 2"),
                    f1 = c(f1(hp$m1, hp$m2), f1(hp$a1, hp$m1), f1(hp$a1, hp$m2), f1(hp$a2, hp$m1), f1(hp$a2, hp$m2)), compared_s = nrow(hp) / 2 / FPS)
for (s in c("m1", "m2", "a1", "a2")) for (cp in c("L", "R")) human[, (paste0(s, "_", cp, "_s")) := sum(hp[cup == cp][[s]]) / FPS]
fwrite(human, file.path(OUT, "socp_frame_human_pair.csv")); print(human, digits = 3)

# ---------------------------------------------------------------- bouts, latency (sessions with the full window)
segs <- function(x, gap, minb) { s <- segment_events(x, FPS, valid = rep(TRUE, length(x)), min_bout_s = minb, max_gap_s = gap)
  c(n = s$n_bouts, mean = s$mean_bout_s, latency = s$latency_s) }
B <- M[full_window == TRUE]
bouts <- rbindlist(lapply(c(0, .2, .5, 1), function(gap) rbindlist(lapply(c(0, .2, .5), function(minb)
  B[, { mm <- segs(m, 0, 0); aa <- segs(a2, gap, minb)
        list(manual_bouts = mm[["n"]], manual_mean_s = mm[["mean"]], manual_latency_s = mm[["latency"]],
             auto_bouts = aa[["n"]], auto_mean_s = aa[["mean"]], auto_latency_s = aa[["latency"]]) },
    by = .(code, phase, cup, role, camera)][, `:=`(merge_gap_s = gap, min_bout_s = minb)]))))
fwrite(bouts, file.path(OUT, "socp_bout_sessions.csv"))
bagree <- bouts[, { ok <- is.finite(manual_mean_s) & is.finite(auto_mean_s); okl <- is.finite(manual_latency_s) & is.finite(auto_latency_s)
  list(cup_sessions = .N, bouts_ratio = sum(auto_bouts) / sum(manual_bouts), bouts_ccc = ccc(auto_bouts, manual_bouts),
       mean_bout_ratio = mean(auto_mean_s[ok]) / mean(manual_mean_s[ok]), mean_bout_ccc = ccc(auto_mean_s[ok], manual_mean_s[ok]),
       latency_ccc = ccc(auto_latency_s[okl], manual_latency_s[okl]), latency_mad_s = mean(abs(auto_latency_s[okl] - manual_latency_s[okl]))) },
  by = .(camera, merge_gap_s, min_bout_s)][order(camera, merge_gap_s, min_bout_s)]
fwrite(bagree, file.path(OUT, "socp_bout_agreement.csv"))
print(bagree[merge_gap_s %in% c(0, .5) & min_bout_s %in% c(0, .2)], digits = 3)
