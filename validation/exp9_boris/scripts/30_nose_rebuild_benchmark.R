# ============================================================================
# 30_nose_rebuild_benchmark.R
#
# How well the landmark repair rebuilds a missing nose, measured where the
# truth is known, with several time constants for the pose of the moment
# (landmark_repair_pose_decay_s; 0 = the typical pose only):
#   hidden stretches   stretches of observed nose (10, 3, 1 and 0.3 s; up to
#                      ten of each per recording, longest placed first, kept
#                      3 s apart) are hidden and rebuilt;
#   real gap edges     the 6 observed frames on either side of every real gap
#                      (nose missing or removed, at least 3 frames) are hidden,
#                      so the rebuilt nose is scored where real gaps begin and
#                      end, under a cup lid included.
# Error = distance between rebuilt and observed nose; near a cup = observed nose
# within 10 cm of a cup point; contact agrees = both within 6 cm of a cup point
# or both beyond (raw noses, before the head length restoration). Every time
# constant is scored on the same frames: those rebuilt under all of them; the
# share each one rebuilt is reported separately. Recordings: every one with at
# least 5 s of rebuilt nose in the landmark-repair release, and 30 others.
#
#   Rscript 30_nose_rebuild_benchmark.R
#   writes results/socp_nose_rebuild_benchmark.csv (per arm, gap length, time constant, near a cup or all)
#   and results/socp_nose_rebuild_benchmark_placed.csv (hidden stretches placed per length)
# ============================================================================

suppressMessages({ library(sp); library(stringr); library(yaml); library(data.table) })
SCRIPTS <- local({ a <- commandArgs(trailingOnly = FALSE); p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/") })
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
REPO <- normalizePath(file.path(PROJ, "..", ".."), winslash = "/")
Sys.setenv(SLEAP_ANALYZER_REPO_ROOT = REPO)
source(file.path(REPO, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))
source(file.path(REPO, "02_SLEAPanalzyer", "core", "assay_config.R"))
OUT <- file.path(PROJ, "results", "socp_nose_rebuild_benchmark.csv")
STAGE <- file.path(Sys.getenv("EXP9_STAGING", "C:/Users/topohl/Documents/exp9_sleap_staging"), "SocP")
REL <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/SLEAPanalyzer_v2/releases/2026-11-05_exp9_all-batches_canonical-landmark-repair"
co <- rbindlist(lapply(list.files(file.path(REL, "assay_summaries", "SocP"), "^combined_output[.]csv$", recursive = TRUE, full.names = TRUE), fread), fill = TRUE)
set.seed(30)
pick <- rbind(co[noseReconstructedS >= 5, .(batch, phase, Code)], co[noseReconstructedS < 5][sample(.N, 30), .(batch, phase, Code)])
DECAY <- c(0, .5, 1, 2); LENGTHS <- c(10, 3, 1, .3); EDGE <- 6L
res <- list(); placed <- list()
for (i in seq_len(nrow(pick))) {
  t <- pick[i]; config <- load_assay_config(file.path(PROJ, "config", if (t$phase == "HAB") "socp_hab_all.yaml" else "socp_all.yaml"), "SocP")
  tr <- read_tracking_csv(file = file.path(STAGE, t$batch, t$phase, paste0(t$Code, "_", t$phase, ".csv")), fps = config$fps)
  tr <- CalibrateTrackingData(tr, method = "area", in.metric = config$arena_width_cm * config$arena_height_cm, points = config$arena_corner_names)
  repair_with <- function(x, decay) repair_landmarks(x, landmarks = config$landmark_repair_landmarks, targets = config$landmark_repair_targets,
    margin = config$landmark_repair_margin_cm, margin_fraction = config$landmark_repair_margin_fraction, min_share = config$landmark_repair_min_share,
    min_landmarks = config$landmark_repair_min_landmarks, max_anchor_error = config$landmark_repair_max_anchor_error_cm,
    corners = config$arena_corner_names, max_outside = config$landmark_repair_outside_cm, pose_decay_s = decay)
  base <- repair_with(tr, 0); st <- as.character(base$data$nose$status)
  ok <- st == "observed" & is.finite(base$data$leftEar$x) & is.finite(base$data$rightEar$x)
  n <- length(ok); fps <- config$fps; margin <- 3 * fps
  cup <- function(p) cbind(stats::median(tr$data[[p]]$x, na.rm = TRUE), stats::median(tr$data[[p]]$y, na.rm = TRUE))
  cl <- cup("socl"); cr <- cup("socr")
  near_d <- function(xy) pmin(sqrt((xy[, 1] - cl[1])^2 + (xy[, 2] - cl[2])^2), sqrt((xy[, 1] - cr[1])^2 + (xy[, 2] - cr[2])^2))
  score <- function(rows, gap_of, arm) {
    truth <- cbind(tr$data$nose$x[rows], tr$data$nose$y[rows]); dt <- near_d(truth)
    masked <- tr; masked$data$nose$x[rows] <- NA; masked$data$nose$y[rows] <- NA
    rbindlist(lapply(DECAY, function(decay) { r <- repair_with(masked, decay); est <- cbind(r$data$nose$x[rows], r$data$nose$y[rows])
      data.table(batch = t$batch, Code = t$Code, phase = t$phase, arm = arm, decay_s = decay, gap_s = gap_of, frame = rows,
                 rebuilt = is.finite(est[, 1]), error_cm = sqrt(rowSums((est - truth)^2)), near_cup = dt <= 10,
                 contact_agrees = (dt <= 6) == (near_d(est) <= 6)) }))
  }
  # hidden stretches: observed nose and ears throughout, 3 s of margin on either side, no overlap, longest first
  taken <- logical(n); H <- list()
  for (len in LENGTHS) {
    L <- round(len * fps); run <- stats::filter(as.numeric(ok), rep(1, L + 2 * margin), sides = 1)
    cand <- which(run == L + 2 * margin) - L - 2 * margin + 1 + margin; got <- 0
    for (s in cand[sample.int(length(cand))]) {
      if (got >= 10) break
      w <- (s - margin):(s + L - 1 + margin)
      if (any(taken[w])) next
      taken[w] <- TRUE; H[[length(H) + 1]] <- data.table(gap_s = len, frame = s:(s + L - 1)); got <- got + 1
    }
    placed[[length(placed) + 1]] <- data.table(batch = t$batch, Code = t$Code, phase = t$phase, gap_s = len, stretches = got)
  }
  if (length(H)) { H <- rbindlist(H); res[[length(res) + 1]] <- score(H$frame, H$gap_s, "hidden stretches") }
  # real gap edges: the EDGE observed frames on either side of each real gap of at least 3 frames, kept apart from other gaps
  r <- rle(st != "observed"); end <- cumsum(r$lengths); start <- end - r$lengths + 1
  gaps <- which(r$values & r$lengths >= 3)
  edges <- unlist(lapply(gaps, function(g) { e <- c((start[g] - EDGE):(start[g] - 1), (end[g] + 1):(end[g] + EDGE)); e <- e[e >= 1 & e <= n]
    if (length(e) == 2 * EDGE && all(ok[e])) e else integer() }))
  edges <- edges[!duplicated(edges)]
  if (length(edges)) res[[length(res) + 1]] <- score(edges, NA_real_, "real gap edges")
  if (i %% 10 == 0) cat(i, "of", nrow(pick), "\n")
}
R <- rbindlist(res)
common <- R[, .(all_rebuilt = all(rebuilt)), by = .(Code, phase, arm, frame)]
R <- common[R, on = .(Code, phase, arm, frame)]
summarise <- function(D, subset) {
  share <- D[, .(rebuilt_share = mean(rebuilt)), by = .(arm, gap_s, decay_s)]
  m <- D[all_rebuilt == TRUE, .(frames = .N, error_median_cm = stats::median(error_cm), error_p95_cm = stats::quantile(error_cm, .95),
                                within_1cm = mean(error_cm <= 1), contact_agrees = mean(contact_agrees)), by = .(arm, gap_s, decay_s)]
  cbind(subset = subset, share[m, on = .(arm, gap_s, decay_s)])
}
summ <- rbind(summarise(R, "all"), summarise(R[near_cup == TRUE], "near a cup"))
summ[, `:=`(animals = uniqueN(R$Code), recordings = uniqueN(R[, paste(Code, phase)]))]
fwrite(summ[order(arm, subset, gap_s, decay_s)], OUT)
fwrite(rbindlist(placed)[, .(recordings = sum(stretches > 0), stretches = sum(stretches)), by = gap_s], sub("[.]csv$", "_placed.csv", OUT))
print(summ[subset == "near a cup"][order(arm, gap_s, decay_s)], digits = 3)
