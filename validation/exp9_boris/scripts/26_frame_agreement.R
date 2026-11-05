# ============================================================================
# 26_frame_agreement.R
#
# The automated social-preference contact against the manual scoring frame by
# frame, and the head length restoration it led to. Batch 1 is the only batch
# scored by hand; for the 10 animals of its second camera the raw BORIS
# project has every event (interaction, rearing, climbing per cup), for the
# other 10 only the transcribed totals.
#
#   Rscript 26_frame_agreement.R <reference_release>
#
# reference_release: the SLEAPanalyzer release the automated values without
# the restoration must reproduce (its assay_summaries/SocP/B1).
#
# Each recording is processed as DLCA_SocP does (calibration, landmark repair,
# interpolation, on-cup detection), and contact is scored twice: with the nose
# as tracked and with the head length restored (restore_head_length()). The
# manual events are put on the video frames and aligned per recording (the
# offset with the best frame agreement within 5 s). OQ752 (A1L9) is left out:
# its phase-2 record describes its phase-1 trial (manual_exclusions.csv of the
# validation folder) and its phase-1 record finds no offset inside the window.
#
# Writes results/socp_frame_*.csv:
#   socp_frame_summary.csv    frame agreement (F1), missed and extra time, totals r / CCC / bias, both noses
#   socp_frame_by_side.csv    manual contact by side of the cup: nose distance and share missed
#   socp_frame_lid_edge.csv   the place along the lid edge: visible head length, scored and unscored
#   socp_frame_sweep.csv      frame agreement and totals against the contact distance, both noses
#   socp_frame_distance.csv   nose distance histogram, frames with and without manual contact
#   socp_frame_positions.csv  nose positions around the cup (0.5 cm grid), by agreement, both noses
#   socp_frame_totals.csv     all 20 animals: automated totals both ways against the manual totals
#   socp_frame_totals_agreement.csv   their r / CCC / bias / mean absolute difference, by measure and camera
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
PAIRS <- file.path(E9, "Analysis/Behavior/Validation/stats/validation_pairs.csv")
OUT <- file.path(PROJ, "results")
config <- load_assay_config(file.path(PROJ, "config", "socp_all.yaml"), "SocP")
stopifnot("the configuration restores the head length" = isTRUE(config$head_length_restore))
qcDecisions <- read_qc_decisions(config$qc_decisions_file)
FPS_VIDEO <- 29.97                       # BORIS times are seconds of the video; tracking rows are its frames
GX <- c("0786" = "A1L9", "0787" = "C7E2", "0788" = "F9L3", "0789" = "K4T2", "0790" = "M1M3", "0791" = "P9X0", "0792" = "R8O3",
        "0793" = "T2H7", "0794" = "V2P0", "0795" = "Y0F3")
LEFT_OUT <- "A1L9"
ALL20 <- sub("_S1[.]csv$", "", list.files(file.path(STAGE, "S1"), "_S1[.]csv$"))
stopifnot(length(ALL20) == 20, all(GX %in% ALL20))
dist <- function(x1, y1, x2, y2) sqrt((x1 - x2)^2 + (y1 - y2)^2)
ccc <- function(x, y) 2 * stats::cov(x, y) / (stats::var(x) + stats::var(y) + (mean(x) - mean(y))^2)

# ---------------------------------------------------------------- per frame, as DLCA_SocP
frames <- rbindlist(lapply(c("S1", "S2"), function(ph) rbindlist(lapply(ALL20, function(code) {
  name <- paste0(code, "_", ph)
  tr <- read_tracking_csv(file = file.path(STAGE, ph, paste0(name, ".csv")), fps = config$fps)
  tr <- CalibrateTrackingData(tr, method = "area", in.metric = config$arena_width_cm * config$arena_height_cm, points = config$arena_corner_names)
  tr <- repair_landmarks(tr, landmarks = config$landmark_repair_landmarks, targets = config$landmark_repair_targets,
    margin = config$landmark_repair_margin_cm, margin_fraction = config$landmark_repair_margin_fraction, min_share = config$landmark_repair_min_share,
    min_landmarks = config$landmark_repair_min_landmarks, max_anchor_error = config$landmark_repair_max_anchor_error_cm,
    corners = config$arena_corner_names, max_outside = config$landmark_repair_outside_cm)
  tr <- interpolate_tracking(tr, landmarks = c("nose", "bodycentre"), max_gap_s = config$max_interpolation_gap_s, likelihood_cutoff = config$likelihood_cutoff)
  size <- apparent_body_size(tr)
  on <- lapply(c(left = "socl", right = "socr"), function(cup) {
    f <- detect_on_object(tr, cup, config$arena_corner_names, footprint_radius = config$on_cup_footprint_radius_cm, beyond = config$on_cup_beyond_cm,
      max_distance = config$on_cup_max_distance_cm, min_apparent_size = config$on_cup_min_apparent_size, merge_gap_s = config$on_cup_merge_gap_s,
      min_episode_s = config$on_cup_min_episode_s, size_window_s = config$on_cup_size_window_s, size = size)
    apply_qc_decisions(f$on_object, qcDecisions, name, cup, config$fps) })
  score <- function(t) compute_socp_metrics(t, "R", config$fps, contact_distance = config$contact_distance_cm,
    body_exclusion_distance = config$body_exclusion_distance_cm, proximity_range = config$proximity_range_cm, threshold_unit = "cm", on_cup = on)$summary
  before <- score(tr); valid <- landmark_validity(tr, c("nose", "bodycentre", "socl", "socr"))
  g <- function(t, p, xy) t$data[[p]][[xy]]
  rs <- restore_head_length(tr); after <- score(rs)
  data.table(code = code, phase = ph, frame = seq_along(g(tr, "nose", "x")), valid = valid, on_l = on$left, on_r = on$right,
             nose_x = g(tr, "nose", "x"), nose_y = g(tr, "nose", "y"), head_x = g(rs, "nose", "x"), head_y = g(rs, "nose", "y"),
             body_x = g(tr, "bodycentre", "x"), body_y = g(tr, "bodycentre", "y"), socl_x = g(tr, "socl", "x"), socl_y = g(tr, "socl", "y"),
             socr_x = g(tr, "socr", "x"), socr_y = g(tr, "socr", "y"), head_ratio = rs$head_length$length / rs$head_length$median,
             left_before = before$contactLeft, right_before = before$contactRight, left_after = after$contactLeft, right_after = after$contactRight)
}))))
rel <- rbindlist(lapply(c("S1", "S2"), function(ph) fread(file.path(REL, "assay_summaries", "SocP", "B1", ph, "combined_output.csv"))))
chk <- unique(frames[, .(file = paste0(code, "_", phase), left_before, right_before)])[rel, on = "file"]
stopifnot("without the restoration the release values are reproduced" = nrow(chk) == 40,
          max(abs(chk$left_before - chk$contactLeft), abs(chk$right_before - chk$contactRight)) < 1e-6)
cat("40 recordings processed; the release values are reproduced without the restoration\n")

# long by cup; contact as compute_socp_metrics scores it
nl <- rbindlist(lapply(c("S1", "S2"), function(ph) fread(sprintf(NOVEL, ph))[, .(code = Code, phase = ph, NovelLoc)]))
frames <- nl[frames, on = .(code, phase)]
cupwise <- function(side) {
  s <- if (side == "L") "socl" else "socr"; o <- if (side == "L") "socr" else "socl"
  frames[, .(code, phase, frame, cup = side, role = fifelse((NovelLoc == "R") == (side == "L"), "mouse", "empty or familiar"),
             valid, on = get(paste0("on_", tolower(side))), head_ratio,
             d = dist(nose_x, nose_y, get(paste0(s, "_x")), get(paste0(s, "_y"))), dh = dist(head_x, head_y, get(paste0(s, "_x")), get(paste0(s, "_y"))),
             b = dist(body_x, body_y, get(paste0(s, "_x")), get(paste0(s, "_y"))),
             rx = (nose_x - get(paste0(s, "_x"))) * sign(get(paste0(o, "_x")) - get(paste0(s, "_x"))), ry = nose_y - get(paste0(s, "_y")),
             hx = (head_x - get(paste0(s, "_x"))) * sign(get(paste0(o, "_x")) - get(paste0(s, "_x"))), hy = head_y - get(paste0(s, "_y")))]
}
L <- rbind(cupwise("L"), cupwise("R"))
L[, ok := valid & !on & b > config$body_exclusion_distance_cm]
L[, `:=`(aut = ok & d <= config$contact_distance_cm, auth = ok & dh <= config$contact_distance_cm)]
L[is.na(aut), aut := FALSE]; L[is.na(auth), auth := FALSE]

# ---------------------------------------------------------------- manual events on the frames
obs <- fromJSON(BORIS, simplifyVector = FALSE)$observations
manual_frames <- function(g, o, n, offset) {
  ob <- obs[[paste0("GX01", g, "_", o)]]; t0 <- ob[["observation time interval"]][[1]]
  t <- vapply(ob$events, function(z) as.numeric(z[[1]]), 0) - t0; b <- vapply(ob$events, function(z) z[[3]], "")
  fr <- seq_len(n) - 1
  state <- function(beh) { k <- which(b == beh); stopifnot(length(k) %% 2 == 0); s <- t[k[c(TRUE, FALSE)]]; e <- t[k[c(FALSE, TRUE)]]; x <- logical(n)
    for (i in seq_along(s)) x[fr >= round((s[i] + offset) * FPS_VIDEO) & fr < round((e[i] + offset) * FPS_VIDEO)] <- TRUE; x }
  list(L = state("l interact") | state("l rear"), R = state("r interact") | state("r rear"), rearL = state("l rear"), rearR = state("r rear"))
}
M <- L[code %in% setdiff(GX, LEFT_OUT)]
setorder(M, code, phase, cup, frame)
lags <- list()
for (g in names(GX)[GX != LEFT_OUT]) for (o in c("stim1", "stim2")) {
  ph <- if (o == "stim1") "S1" else "S2"; k <- M$code == GX[[g]] & M$phase == ph; n <- sum(k) / 2
  aL <- M$aut[k & M$cup == "L"]; aR <- M$aut[k & M$cup == "R"]
  off <- seq(-5, 5, by = 1 / 30)
  agree <- vapply(off, function(z) { m <- manual_frames(g, o, n, z); sum(m$L == aL) + sum(m$R == aR) }, 0)
  best <- off[which.max(agree)]; stopifnot("the alignment is inside the search window" = abs(best) < 4.9)
  m <- manual_frames(g, o, n, best)
  M[k & cup == "L", `:=`(man = m$L, rear = m$rearL)]; M[k & cup == "R", `:=`(man = m$R, rear = m$rearR)]
  lags[[length(lags) + 1]] <- data.table(code = GX[[g]], phase = ph, offset_s = best)
}
lags <- rbindlist(lags)
cat(sprintf("manual events aligned: %d recordings, offsets %.2f to %.2f s (median %.2f)\n", nrow(lags), min(lags$offset_s), max(lags$offset_s), median(lags$offset_s)))

# ---------------------------------------------------------------- frame agreement
agree_tab <- function(pred, by = NULL) {
  r <- M[, .(tp = sum(man & get(pred)), fn = sum(man & !get(pred)), fp = sum(!man & get(pred)), man_s = sum(man) / 30, aut_s = sum(get(pred)) / 30),
         by = c("code", "phase", "cup", "role")]
  r[, .(cup_recordings = .N, f1 = 2 * sum(tp) / (2 * sum(tp) + sum(fn) + sum(fp)), manual_s = sum(man_s), both_s = sum(tp) / 30,
        missed_s = sum(fn) / 30, extra_s = sum(fp) / 30, r = stats::cor(aut_s, man_s), ccc = ccc(aut_s, man_s),
        bias_s = mean(aut_s - man_s), loa_s = 1.96 * stats::sd(aut_s - man_s)), by = by]
}
summ <- rbind(cbind(nose = "as tracked", agree_tab("aut")), cbind(nose = "head length restored", agree_tab("auth")))
fwrite(summ, file.path(OUT, "socp_frame_summary.csv")); print(summ, digits = 3)

# by side of the cup (angle of the nose around the cup point; 0 = toward the other cup)
M[, angle := abs(((atan2(ry, rx) * 180 / pi) + 180) %% 360 - 180)]
M[, side := cut(angle, c(0, 45, 90, 135, 180), include.lowest = TRUE, labels = c("0-45", "45-90", "90-135", "135-180"))]
side <- M[man & valid, .(manual_s = .N / 30, median_nose_cm = stats::median(d, na.rm = TRUE), missed_as_tracked = mean(!aut), missed_head_restored = mean(!auth)),
          by = .(role, side)][order(role, side)]
fwrite(side, file.path(OUT, "socp_frame_by_side.csv")); print(side, digits = 3)

# along the lid edge: beyond the contact distance, on the wall side, behind the cup point
M[, lid_edge := valid & d > 6.3 & d < 11 & abs(ry) > 3.8 & rx < -2.5]
lid <- M[lid_edge == TRUE & is.finite(head_ratio), .(time_s = .N / 30, median_head_ratio = stats::median(head_ratio), share_head_below_75 = mean(head_ratio < .75),
                                                     contact_head_restored = mean(auth)), by = .(manual_contact = man)]
fwrite(lid, file.path(OUT, "socp_frame_lid_edge.csv")); print(lid, digits = 3)

# contact distance sweep, both noses
sweep <- rbindlist(lapply(seq(4, 9, by = .25), function(cd) rbindlist(lapply(c(`as tracked` = "d", `head length restored` = "dh"), function(v) {
  M[, z := ok & get(v) <= cd]; M[is.na(z), z := FALSE]
  cbind(contact_distance_cm = cd, agree_tab("z")) }), idcol = "nose")))
M[, z := NULL]
fwrite(sweep, file.path(OUT, "socp_frame_sweep.csv"))

# nose distance histogram and positions around the cup
br <- seq(0, 20, by = .25)
hist <- rbindlist(lapply(c(`as tracked` = "d", `head length restored` = "dh"), function(v)
  M[valid & get(v) < 20, .(time_s = .N / 30), by = .(manual_contact = man, bin_cm = br[findInterval(get(v), br)] + .125)]), idcol = "nose")
fwrite(hist[order(nose, manual_contact, bin_cm)], file.path(OUT, "socp_frame_distance.csv"))
st <- function(m, a) fifelse(m & a, "both", fifelse(m, "manual only", fifelse(a, "automated only", "neither")))
pos <- rbind(M[valid & d < 14, .(nose = "as tracked", x = round(rx * 2) / 2, y = round(abs(ry) * 2) / 2, status = st(man, aut))],
             M[valid & dh < 14, .(nose = "head length restored", x = round(hx * 2) / 2, y = round(abs(hy) * 2) / 2, status = st(man, auth))])
fwrite(pos[status != "neither", .(time_s = .N / 30), by = .(nose, status, x, y)], file.path(OUT, "socp_frame_positions.csv"))

# ---------------------------------------------------------------- totals, all 20 animals
pairs <- fread(PAIRS)[assay == "SocP" & unit == "s" & is.finite(boris)]
pairs[, `:=`(phase = substr(metric, 1, 2), role = fifelse(grepl("novel", metric), "mouse", "empty or familiar"))]
tot <- L[, .(as_tracked = sum(aut) / 30, head_restored = sum(auth) / 30), by = .(code, phase, role)]
tot <- pairs[, .(code = Code, Animal, phase, role, metric, manual = boris)][tot, on = .(code, phase, role), nomatch = 0]
tot[, camera := fifelse(code %in% GX, "camera 2 (raw events)", "camera 1 (totals only)")]
fwrite(tot, file.path(OUT, "socp_frame_totals.csv"))
ta <- function(x) x[, .(n = .N, r_as_tracked = stats::cor(as_tracked, manual), r_head_restored = stats::cor(head_restored, manual),
                        ccc_as_tracked = ccc(as_tracked, manual), ccc_head_restored = ccc(head_restored, manual),
                        bias_as_tracked = mean(as_tracked - manual), bias_head_restored = mean(head_restored - manual),
                        mad_as_tracked = mean(abs(as_tracked - manual)), mad_head_restored = mean(abs(head_restored - manual)))]
tagree <- rbind(tot[, ta(.SD), by = .(metric)][, camera := "both"], tot[, ta(.SD), by = .(metric, camera)])
fwrite(tagree, file.path(OUT, "socp_frame_totals_agreement.csv")); print(tagree, digits = 3)
