# ============================================================================
# 24_qc_report.R
#
# The tracking QC report of the Exp9 social preference test, as tables (one
# formatted workbook, the same sheets as CSV) and figures (SVG, PNG preview):
# what the frame-level checks found, what was repaired or left out, and what
# that changed. Report only: it reads the SocP outputs and changes nothing.
#
#   Rscript 24_qc_report.R <socp_output_root> <reference_release> <out_dir>
#
#   socp_output_root   folder with <batch>/<phase>/combined_output.csv and the
#                      QC tables of DLCA_SocP (on_cup_episodes.csv,
#                      landmark_repairs.csv, landmark_repair_episodes.csv,
#                      landmark_repair_models.csv)
#   reference_release  the SLEAPanalyzer release before the repair, for the
#                      before/after comparison (its assay_summaries/SocP)
#   out_dir            written: qc_report.xlsx, tables/*.csv, figures/*.svg|png
#
# Also read: results/plausibility_summary.csv (23_plausibility_screen.R), the
# reviewed decisions (config/qc_decisions_socp.csv), the release animal table
# for sex and group, and the batch 1 manual scores
# (Analysis/Behavior/Validation/stats/validation_pairs.csv) when present.
# ============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(openxlsx); library(yaml) })
args <- commandArgs(TRUE); if (length(args) != 3) stop("usage: Rscript 24_qc_report.R <socp_output_root> <reference_release> <out_dir>")
ROOT <- args[1]; REF <- args[2]; OUT <- args[3]
SCRIPTS <- local({ a <- commandArgs(trailingOnly = FALSE); p <- sub("^--file=", "", a[grep("^--file=", a)]); normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/") })
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
A <- Sys.getenv("E9_ANALYSIS_ROOT", "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis")
PAL <- yaml::read_yaml(file.path(A, "Exp9_manuscript", "config", "manuscript_palette.yml"))
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)
BATCHES <- paste0("B", 1:6); PHASES <- c("HAB", "S1", "S2")
PHASE_LABEL <- c(HAB = "Habituation (empty cups)", S1 = "Phase 1 (mouse vs empty cup)", S2 = "Phase 2 (new vs familiar mouse)")
rd_all <- function(root, name) rbindlist(lapply(BATCHES, function(b) rbindlist(lapply(PHASES, function(p) {
  f <- file.path(root, b, p, name); if (!file.exists(f) || file.size(f) < 5) return(NULL); x <- fread(f); if (!nrow(x)) return(NULL); x[, `:=`(batch = b, phase = p)] }), fill = TRUE)), fill = TRUE)

# ---------------------------------------------------------------- inputs
now <- rd_all(ROOT, "combined_output.csv"); ref <- rd_all(file.path(REF, "assay_summaries", "SocP"), "combined_output.csv")
reps <- rd_all(ROOT, "landmark_repairs.csv"); repe <- rd_all(ROOT, "landmark_repair_episodes.csv"); mods <- rd_all(ROOT, "landmark_repair_models.csv")
cups <- rd_all(ROOT, "on_cup_episodes.csv")
hand <- fread(file.path(REF, "data", "exp9_sleap_all_animals.tsv"))
meta <- hand[, .(Code, ID, Sex, Group = Phenotype)]
dec <- fread(file.path(PROJ, "config", "qc_decisions_socp.csv"))
screen_f <- file.path(PROJ, "results", "plausibility_summary.csv"); screen <- if (file.exists(screen_f)) fread(screen_f) else NULL
stopifnot("every recording also in the reference release" = all(paste(now$batch, now$phase, now$Code) %in% paste(ref$batch, ref$phase, ref$Code)),
          "every SocP output comes from a run with the landmark repair" = all(c("landmarksFlaggedS", "noseReconstructedS") %in% names(now)) &&
            !anyNA(now$landmarksFlaggedS) && !anyNA(now$noseReconstructedS),
          "a repair table for every recording" = all(paste(now$batch, now$phase, now$Code) %in% paste(reps$batch, reps$phase, reps$Code)))
commits <- unique(unlist(lapply(BATCHES, function(b) lapply(PHASES, function(p) { f <- file.path(ROOT, b, p, "run_manifest.yaml")
  if (file.exists(f)) yaml::read_yaml(f)$code$commit }))))
cat("SocP outputs:", nrow(now), "recordings; SLEAPanalyzer commit(s):", paste(substr(commits, 1, 7), collapse = ", "), "\n")

# ---------------------------------------------------------------- recordings
k <- match(paste(now$batch, now$phase, now$Code), paste(ref$batch, ref$phase, ref$Code))
j <- match(now$Code, meta$Code); stopifnot("every recording has an animal in the release table" = !anyNA(j))
rec <- data.table(batch = now$batch, phase = now$phase, Code = now$Code, ID = meta$ID[j], Sex = meta$Sex[j], Group = meta$Group[j])
rec[, `:=`(valid_fraction = now$valid_fraction, qc_pass = now$qcPass,
           landmarks_flagged_s = now$landmarksFlaggedS, skeleton_collapse_s = now$skeletonCollapseS,
           nose_rebuilt_s = now$noseReconstructedS, bodycentre_rebuilt_s = now$bodycentreReconstructedS,
           nose_anchor_error_p95_cm = now$noseAnchorErrorP95,
           on_cup_s = now$onCupLeft + now$onCupRight, on_cup_contact_removed_s = now$contactRemovedLeft + now$contactRemovedRight,
           novel_before_s = ref$contactNovel[k], novel_after_s = now$contactNovel,
           familiar_before_s = ref$contactFamiliar[k], familiar_after_s = now$contactFamiliar,
           distance_before_cm = ref$distance[k], distance_after_cm = now$distance)]
rec[, `:=`(novel_change_s = novel_after_s - novel_before_s, familiar_change_s = familiar_after_s - familiar_before_s,
           preference_before = (novel_before_s - familiar_before_s) / (novel_before_s + familiar_before_s),
           preference_after = (novel_after_s - familiar_after_s) / (novel_after_s + familiar_after_s))]
rec[, preference_change := preference_after - preference_before]
split <- if (nrow(repe)) repe[check == "skeleton_split", .(skeleton_split_s = sum(duration_s)), by = .(batch, phase, Code)] else data.table(batch = character(), phase = character(), Code = character(), skeleton_split_s = numeric())
rec <- split[rec, on = c("batch", "phase", "Code")]; rec[is.na(skeleton_split_s), skeleton_split_s := 0]
setcolorder(rec, c(setdiff(names(rec)[seq_len(match("skeleton_collapse_s", names(rec)))], "skeleton_split_s"), "skeleton_split_s"))
setorder(rec, batch, phase, Code)
reps <- meta[reps, on = "Code"]; repe <- meta[repe, on = "Code"]; mods <- meta[mods, on = "Code"]; cups <- meta[cups, on = "Code"]
if (nrow(cups)) {
  cups[, reviewed := ""]
  for (i in seq_len(nrow(dec))) cups[paste0(Code, "_", phase) == dec$file[i] & object == dec$target[i] & start_s <= dec$end_s[i] & start_s + duration_s >= dec$start_s[i],
                                     reviewed := paste0(dec$decision[i], ": ", dec$reason[i])]
}

# ---------------------------------------------------------------- manual check (batch 1)
man <- NULL
vp_f <- file.path(A, "Behavior", "Validation", "stats", "validation_pairs.csv")
if (file.exists(vp_f)) {
  vp <- fread(vp_f)[assay == "SocP" & metric %in% c("S1 novel", "S1 familiar", "S2 novel", "S2 familiar")]
  vp[, `:=`(phase = substr(metric, 1, 2), side = sub("^S[12] ", "", metric))]
  get_val <- function(d, ph, sd, code) { x <- d[phase == ph & Code == code]; if (!nrow(x)) NA_real_ else if (sd == "novel") x$contactNovel else x$contactFamiliar }
  vp[, before := mapply(get_val, list(ref), phase, side, Code)]; vp[, after := mapply(get_val, list(now), phase, side, Code)]
  man <- vp[, .(ID = Animal, Code, phase, side, manual = boris, manual_as_transcribed = if ("boris_as_transcribed" %in% names(vp)) boris_as_transcribed else boris,
                excluded = if ("excluded" %in% names(vp)) excluded else "", sleap_before = before, sleap_after = after)]
  ccc <- function(x, y) { n <- length(x); vx <- var(x) * (n - 1) / n; vy <- var(y) * (n - 1) / n; cxy <- cov(x, y) * (n - 1) / n; 2 * cxy / (vx + vy + (mean(x) - mean(y))^2) }
  agree <- man[is.finite(manual), .(n = .N, r_before = cor(sleap_before, manual), r_after = cor(sleap_after, manual), ccc_before = ccc(sleap_before, manual),
                                    ccc_after = ccc(sleap_after, manual), bias_before = mean(sleap_before - manual), bias_after = mean(sleap_after - manual)), by = .(phase, side)]
}

# ---------------------------------------------------------------- tables
r3 <- function(d) { d <- copy(d); for (cn in names(d)) if (is.numeric(d[[cn]]) && !is.integer(d[[cn]])) set(d, j = cn, value = round(d[[cn]], 3)); d }
TAB <- list(
  Recordings = rec,
  `Landmark repairs` = reps[, c("ID", "Code", "Sex", "Group", "batch", "phase", "landmark", "outside_s", "flagged_s", "collapse_s", "reconstructed_s", "still_missing_s"), with = FALSE],
  `Repair episodes` = repe[, intersect(c("ID", "Code", "Sex", "Group", "batch", "phase", "check", "target", "start_s", "duration_s", "reconstructed_share", "start_frame", "end_frame"), names(repe)), with = FALSE],
  `On-cup episodes` = if (nrow(cups)) cups[, intersect(c("ID", "Code", "Sex", "Group", "batch", "phase", "object", "start_s", "duration_s", "decision", "reviewed", "median_apparent_size", "min_distance", "max_beyond", "start_frame", "end_frame"), names(cups)), with = FALSE] else data.table(),
  `Anchor models` = mods[, intersect(c("ID", "Code", "batch", "phase", "target", "rank", "anchor_a", "anchor_b", "frames", "error_median", "error_p95"), names(mods)), with = FALSE],
  `Reviewed decisions` = dec)
if (!is.null(man)) { TAB$`Manual check` <- man; TAB$`Manual agreement` <- agree }
if (!is.null(screen)) TAB$`Screen, all tests` <- screen
for (nm in names(TAB)) fwrite(r3(TAB[[nm]]), file.path(OUT, "tables", paste0(gsub("[^A-Za-z0-9]+", "_", tolower(nm)), ".csv")))

n_rep <- rec[nose_rebuilt_s + bodycentre_rebuilt_s > 0, .N]; n_col <- rec[skeleton_collapse_s > 0, .N]; n_cup <- rec[on_cup_s > 0, .N]
readme <- c(
  "Tracking QC report, social preference test (SLEAP, all animals, habituation and both test phases). Report only; the values analysed are in the SLEAPanalyzer release.",
  sprintf("Built %s from %s (reference before the repair: %s).", format(Sys.Date()), ROOT, basename(REF)),
  "",
  "What the checks do (SLEAPanalyzer docs/assay_definitions.md):",
  "- On top of a cup: body centre where a floor-level mouse cannot be and the mouse imaged larger; that time is neither contact nor proximity.",
  "- Landmark repair: a body landmark more than 15 cm outside the arena, or one that disagrees with most of the skeleton (distance to the others far above the recording's typical distance), is removed; when about half the skeleton sits elsewhere, the frame keeps the consistent majority (skeleton split); nose and body centre are rebuilt from two rigid neighbours of the same frame (their typical relative position); frames without a consistent majority (the mouse hidden, hands in the arena) stay untracked (skeleton collapse).",
  "- Reviewed decisions override the automatic rule (sheet Reviewed decisions).",
  "",
  sprintf("Recordings: %d. With a rebuilt nose or body centre: %d. With skeleton collapse: %d. With time on top of a cup: %d.", nrow(rec), n_rep, n_col, n_cup),
  sprintf("Rebuilt nose in total %.0f s, body centre %.0f s; skeleton collapse %.0f s; on top of a cup %.0f s (contact removed %.0f s).",
          sum(rec$nose_rebuilt_s), sum(rec$bodycentre_rebuilt_s), sum(rec$skeleton_collapse_s), sum(rec$on_cup_s), sum(rec$on_cup_contact_removed_s)),
  "",
  "Sheets:",
  "- Recordings: one row per recording: seconds flagged, rebuilt, collapsed, on a cup; the error of the nose rebuild (95th percentile, cm); contact and preference before and after the repair.",
  "- Landmark repairs: per recording and landmark. Repair episodes: every stretch of a landmark outside the arena or inconsistent with the skeleton, and of a split or collapsed skeleton. On-cup episodes: every candidate with the rule's decision and any review.",
  "- Anchor models: the three best neighbour pairs per rebuilt landmark and recording, with their error on clean frames.",
  "- Manual check / Manual agreement (batch 1): manual scores against the automated ones before and after the repair.",
  "- Screen, all tests: every check on every recording of the four tests (23_plausibility_screen.R; report only; NOR, OFT and EPM are not repaired).",
  "",
  "Figures (figures/): repairs_by_recording, contact_change, manual_agreement, landmark_heatmap, anchor_error, screen_all_tests. Females above males throughout.")
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#E6E4DF", border = "bottom", wrapText = TRUE, valign = "top")
addWorksheet(wb, "Read me"); writeData(wb, "Read me", data.frame(`Read me` = readme, check.names = FALSE)); setColWidths(wb, "Read me", 1, 160)
addStyle(wb, "Read me", hs, rows = 1, cols = 1)
for (nm in names(TAB)) { sh <- substr(nm, 1, 31); d <- r3(TAB[[nm]]); addWorksheet(wb, sh); writeData(wb, sh, d, headerStyle = hs)
  freezePane(wb, sh, firstRow = TRUE); setColWidths(wb, sh, seq_len(ncol(d)), "auto")
  num <- which(vapply(d, is.numeric, TRUE))
  for (cn in intersect(c("nose_rebuilt_s", "bodycentre_rebuilt_s", "skeleton_collapse_s", "on_cup_s", "flagged_s", "reconstructed_s", "collapse_s", "duration_s"), names(d)))
    conditionalFormatting(wb, sh, cols = match(cn, names(d)), rows = 2:(nrow(d) + 1), style = c("#FFFFFF", "#F2CA4E"), type = "colourScale")
  for (cn in intersect(c("novel_change_s", "familiar_change_s", "preference_change"), names(d)))
    conditionalFormatting(wb, sh, cols = match(cn, names(d)), rows = 2:(nrow(d) + 1), style = c("#6679D9", "#FFFFFF", "#F2CA4E"), type = "colourScale")
}
saveWorkbook(wb, file.path(OUT, "qc_report.xlsx"), overwrite = TRUE)

# ---------------------------------------------------------------- figures
T <- PAL$typography; GC <- unlist(PAL$group)
th <- theme_classic(base_family = T$family, base_size = T$axis_text_pt) +
  theme(axis.line = element_line(linewidth = PAL$line$axis_pt / 2.13), axis.ticks = element_line(linewidth = PAL$line$axis_pt / 2.13),
        axis.text = element_text(size = T$axis_text_pt), axis.title = element_text(size = T$axis_title_pt), strip.background = element_blank(),
        strip.text = element_text(size = T$axis_title_pt, hjust = 0), legend.text = element_text(size = T$legend_text_pt),
        legend.title = element_text(size = T$legend_title_pt), legend.key.size = unit(3, "mm"), panel.spacing = unit(2, "mm"))
save_fig <- function(p, name, w, h) { ggsave(file.path(OUT, "figures", paste0(name, ".svg")), p, width = w, height = h, units = "mm")
  ggsave(file.path(OUT, "figures", paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white") }
rec[, sex := factor(Sex, c("Female", "Male"))]; rec[, phase_label := factor(PHASE_LABEL[phase], PHASE_LABEL)]; rec[, Group := factor(Group, names(GC))]

long <- melt(rec, id.vars = c("Code", "ID", "sex", "phase_label", "Group", "batch"),
             measure.vars = c("nose_rebuilt_s", "bodycentre_rebuilt_s", "skeleton_split_s", "skeleton_collapse_s", "on_cup_s"), variable.name = "what", value.name = "s")
long[, what := factor(what, c("nose_rebuilt_s", "bodycentre_rebuilt_s", "skeleton_split_s", "skeleton_collapse_s", "on_cup_s"),
                      c("Nose rebuilt", "Body centre rebuilt", "Skeleton split (majority kept)", "Skeleton collapsed (untracked)", "On top of a cup"))]
p1 <- ggplot(long, aes(s, what, colour = Group)) + geom_point(shape = 16, size = 0.9, alpha = 0.8, position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  scale_colour_manual(values = GC, drop = FALSE) + facet_grid(sex ~ phase_label) + scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100)) +
  labs(x = "Seconds per recording", y = NULL, colour = NULL) + th
save_fig(p1, "repairs_by_recording", 183, 78)

cl <- melt(rec, id.vars = c("Code", "sex", "phase_label", "Group"), measure.vars = list(c("novel_before_s", "familiar_before_s"), c("novel_after_s", "familiar_after_s")),
           variable.name = "side", value.name = c("before", "after"))
# novel / familiar side: in habituation the side of the stranger to come, in phase 1 the mouse and the empty cup
PANELS <- c("Habituation, stranger side", "Habituation, other side", "Phase 1, mouse", "Phase 1, empty cup", "Phase 2, new mouse", "Phase 2, familiar mouse")
cl[, panel := factor(PANELS[(match(as.character(phase_label), PHASE_LABEL) - 1) * 2 + as.integer(side)], PANELS)]
p2 <- ggplot(cl, aes(before, after, colour = Group)) + geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = PAL$line$reference_pt / 2.13, colour = "grey55") +
  geom_point(shape = 16, size = 0.9, alpha = 0.8) + scale_colour_manual(values = GC, drop = FALSE) + facet_grid(sex ~ panel) + coord_equal() +
  labs(x = "Cup contact before the repair (s)", y = "After (s)", colour = NULL) + th
save_fig(p2, "contact_change", 183, 75)

if (!is.null(man)) {
  mm <- man[is.finite(manual)]; mm[, panel := factor(paste(phase, side), c("S1 novel", "S1 familiar", "S2 novel", "S2 familiar"),
                                                      c("Phase 1, mouse", "Phase 1, empty cup", "Phase 2, new mouse", "Phase 2, familiar mouse"))]
  ml <- rbind(mm[, .(panel, manual, sleap = sleap_before, when = "before the repair")], mm[, .(panel, manual, sleap = sleap_after, when = "after the repair")])
  ml[, when := factor(when, c("before the repair", "after the repair"))]
  p3 <- ggplot(mm) + geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = PAL$line$reference_pt / 2.13, colour = "grey55") +
    geom_segment(aes(manual, sleap_before, xend = manual, yend = sleap_after), colour = "grey60", linewidth = 0.25) +
    geom_point(data = ml, aes(manual, sleap, colour = when), shape = 16, size = 0.9) +
    scale_colour_manual(values = c(`before the repair` = "#B9B9B4", `after the repair` = "#1F3D52"), name = NULL) +
    facet_wrap(~panel, nrow = 1, scales = "free") + labs(x = "Manual (s)", y = "Automated (s)") + th + theme(aspect.ratio = 1)
  save_fig(p3, "manual_agreement", 183, 55)
}

hm <- reps[, .(s = mean(outside_s + flagged_s + collapse_s)), by = .(landmark, Sex, phase)]
hm[, landmark := factor(landmark, rev(unique(reps$landmark)))]; hm[, col := factor(paste(Sex, phase), paste(rep(c("Female", "Male"), each = 3), PHASES))]
p4 <- ggplot(hm, aes(col, landmark, fill = s)) + geom_tile(width = PAL$heatmap$tile, height = PAL$heatmap$tile) +
  geom_text(aes(label = sprintf("%.1f", s)), size = T$annotation_pt / .pt, family = T$family) +
  scale_fill_gradient(low = "#FFFFFF", high = "#F2CA4E", name = "s removed\nper recording") + coord_equal() + labs(x = NULL, y = NULL) + th +
  theme(axis.line = element_blank(), axis.ticks = element_blank())
save_fig(p4, "landmark_heatmap", 120, 80)

ae <- mods[target %in% c("nose", "bodycentre") & rank == 1]
p5 <- ggplot(ae, aes(error_p95, target)) + geom_vline(xintercept = 2, linetype = "22", linewidth = PAL$line$reference_pt / 2.13, colour = "grey55") +
  geom_point(shape = 16, size = 0.9, alpha = 0.6, colour = "#1F3D52", position = position_jitter(height = 0.2, width = 0, seed = 2)) +
  labs(x = "Error of the rebuild on clean frames, 95th percentile (cm); dashed: limit", y = NULL) + th
save_fig(p5, "anchor_error", 120, 40)

if (!is.null(screen)) {
  sl <- melt(screen[failed == FALSE], id.vars = c("assay", "phase", "batch", "Code"), measure.vars = c("jump_s", "skeleton_s", "outside_arena_s", "raised_s", "on_object_s"),
             variable.name = "check", value.name = "s")
  sl[, check := factor(check, c("jump_s", "skeleton_s", "outside_arena_s", "raised_s", "on_object_s"), c("Jumps", "Skeleton distorted", "Outside the arena", "Raised", "On a cup or object"))]
  sl[, test := factor(paste(assay, phase), c("SocP HAB", "SocP S1", "SocP S2", "NOR HAB", "NOR NOV", "OFT OFT", "EPM EPM"),
                      c("Social, habituation", "Social, phase 1", "Social, phase 2", "Novel object, habituation", "Novel object, test", "Open field", "Plus maze"))]
  p6 <- ggplot(sl, aes(s, test)) + geom_point(shape = 16, size = 0.7, alpha = 0.6, colour = "#1F3D52", position = position_jitter(height = 0.2, width = 0, seed = 3)) +
    facet_wrap(~check, nrow = 1) + scale_x_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100)) + labs(x = "Seconds per recording", y = NULL) + th
  save_fig(p6, "screen_all_tests", 183, 55)
}
cat("wrote", OUT, ":", length(TAB) + 1, "sheets,", length(list.files(file.path(OUT, "figures"), "[.]svg$")), "figures\n")
