# ============================================================================
# 34_hidden_review.R
#
# Review material for the hidden-contact rule of the Exp9 social preference
# test: every stretch the rule counts as hidden under a cup box (decision
# hidden) and every stretch a review cleared (not_hidden), from the
# hidden_episodes.csv of each batch and phase of a SocP output tree. Writes
# the episode and video lists for tools/qc_review.py, which renders a still
# strip and a short clip of each stretch with its cup drawn as a cross. Renders
# nothing itself.
#
#   Rscript 34_hidden_review.R <socp_output_root> <work_dir> [video_qc_manifest.csv]
#   python -I tools/qc_review.py <work_dir>/episodes.csv <work_dir>/videos.csv <out_dir>
#
#   socp_output_root  <root>/<batch>/<phase>/hidden_episodes.csv, as
#                     DLCA_SocP writes it with hidden_contact on (output_dir of
#                     config/socp_all.yaml and socp_hab_all.yaml)
#   work_dir          written: episodes.csv (episode_id, tracking_file,
#                     start_frame, end_frame, label, target, group, note,
#                     video) and videos.csv (tracking_file, video_path)
#   video_qc_manifest the staged tracking file and the clip of each recording
#                     (32_video_qc_manifest.R; default
#                     results/video_qc_manifest.csv, not kept in git)
#
# hidden_episodes.csv stays on the clip's clock after a leave-out, so its
# start_frame and end_frame are rows of the staged tracking file and frames of
# the clip, as qc_review.py expects. Checks, each stopping the script with the
# recordings named: every batch/phase folder with a combined_output.csv has a
# hidden_episodes.csv (a run without hidden_contact is no review base), and
# every listed recording has one SocP row in the manifest whose tracking file
# is the recording's and exists, and a clip. Runtime: seconds.
# ============================================================================

suppressMessages(library(data.table))
SCRIPTS <- local({ a <- commandArgs(trailingOnly = FALSE); p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/") })
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
REPO <- normalizePath(file.path(PROJ, "..", ".."), winslash = "/")
args <- commandArgs(TRUE)
if (!length(args) %in% 2:3) stop("usage: Rscript 34_hidden_review.R <socp_output_root> <work_dir> [video_qc_manifest.csv]")
ROOT <- args[1]; WORK <- normalizePath(args[2], winslash = "/", mustWork = FALSE)
MANIFEST <- if (length(args) == 3) args[3] else file.path(PROJ, "results", "video_qc_manifest.csv")
if (!dir.exists(ROOT)) stop("no SocP output tree at ", ROOT)
if (!file.exists(MANIFEST)) stop("no video manifest at ", MANIFEST, ": run 32_video_qc_manifest.R first")

# The hidden stretches of every batch and phase.
folders <- CJ(batch = grep("^B[0-9]+$", basename(list.dirs(ROOT, recursive = FALSE)), value = TRUE),
              phase = c("HAB", "S1", "S2"))
folders[, dir := file.path(ROOT, batch, phase)]
folders <- folders[file.exists(file.path(dir, "combined_output.csv"))]
if (!nrow(folders)) stop("no SocP output (combined_output.csv) under ", ROOT)
lacking <- folders[!file.exists(file.path(dir, "hidden_episodes.csv"))]
if (nrow(lacking)) stop("no hidden_episodes.csv (scored without hidden_contact?) in: ",
                        paste(lacking$batch, lacking$phase, collapse = ", "))
episodes <- rbindlist(lapply(file.path(folders$dir, "hidden_episodes.csv"), function(f)
  if (file.size(f) > 5) fread(f, colClasses = list(character = c("file", "Code", "batch", "phase")))), fill = TRUE)
if (!nrow(episodes)) stop("no hidden stretch in any hidden_episodes.csv under ", ROOT)
cat("hidden_episodes.csv:", nrow(folders), "batch/phase folders,", nrow(episodes), "body-fill stretches;",
    "decisions:", paste(names(table(episodes$decision)), table(episodes$decision), collapse = ", "), "\n")
episodes <- episodes[decision %in% c("hidden", "not_hidden")]
setorder(episodes, batch, phase, Code, start_frame)

# The staged tracking file and the clip of each recording.
manifest <- fread(MANIFEST, colClasses = "character")[assay == "SocP", .(batch, phase, Code, tracking, video)]
manifest[, tracking := chartr("\\", "/", tracking)]
recordings <- unique(episodes[, .(batch, phase, Code, file)])
found <- manifest[recordings, on = c("batch", "phase", "Code")]
rows <- found[, .N, by = .(batch, phase, Code)]
bad <- found[!is.na(tracking) & (basename(tracking) != paste0(file, ".csv") | !file.exists(tracking)) |
               is.na(video) | !nzchar(video)]
bad <- unique(rbind(bad[, .(batch, phase, Code)], rows[N != 1, .(batch, phase, Code)]))
if (nrow(bad)) stop("no single staged tracking file and clip in the manifest for: ",
                    paste(bad$batch, bad$phase, bad$Code, collapse = ", "))
episodes <- found[, .(batch, phase, Code, tracking, video)][episodes, on = c("batch", "phase", "Code")]

# The lists for qc_review.py.
episodes[, episode_id := sprintf("hid_%s_%s_%02d", batch, file, seq_len(.N)), by = .(batch, file)]
episodes[, group := ifelse(decision == "hidden", "Hidden: counted as contact with the cup",
                           "Not hidden (reviewed): not counted")]
out <- episodes[, .(
  episode_id, tracking_file = tracking, start_frame, end_frame,
  label = sprintf("%s %s %.1f s", decision, target, duration_s), target, group,
  note = sprintf("%.2f-%.2f s; body %.1f-%.1f cm from the cup; nose lost in %.0f%%", start_s, end_s,
                 body_min_distance, body_max_distance, 100 * nose_lost_share),
  video)]
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
fwrite(out, file.path(WORK, "episodes.csv"))
fwrite(unique(out[, .(tracking_file, video_path = video)]), file.path(WORK, "videos.csv"))

summary <- episodes[, .(stretches = .N, seconds = round(sum(duration_s), 1), recordings = uniqueN(file)), by = .(decision, batch)]
print(summary, row.names = FALSE, class = FALSE)
cat("\nwrote", file.path(WORK, "episodes.csv"), "and videos.csv:", nrow(out), "stretches in",
    uniqueN(out$tracking_file), "recordings. Render the stills and clips with\n\n")
cat(sprintf("  python -I %s %s %s %s\n", file.path(REPO, "tools", "qc_review.py"), file.path(WORK, "episodes.csv"),
            file.path(WORK, "videos.csv"), file.path(WORK, "review")))
