# ============================================================================
# 32_video_qc_manifest.R
#
# The manifest for the video QC (tools/video_qc): one row per staged Exp9
# recording, the same recordings 31_light_checks.R checks, with the clip it
# was tracked from and what the checks need to know about its test. The tool
# compares each tracking file with its clip (frame count and length, fixed
# objects, the mouse frame by frame, the clip within its raw recording);
# 31_light_checks.R joins its flags as "video:<check>:<flag>". Report only.
#
#   video      the file in <batch>/<test folder>/cropped under VIDEOS (batch 2
#              novel object: NOR_B2_cropped; subfolders such as old and temp
#              are not searched) whose name, split at "_" (or any character
#              other than a letter or digit), holds the animal
#              code and the recording's trial (HAB, S1, S2, NOV; none in the
#              open field and plus maze). Batches name clips differently
#              (A1L9_GX0786_HAB, C2P4_GX010797, E9_B1_EPM_GX010017_U5P3,
#              E9_SIS_B2_A0S5_HAB, E9_SIS_B2_A0S5, E9_B4_A0I4_OFT); a part
#              that is only a number marks a re-crop (B4Y7_NOV_1), and a novel
#              object clip with no trial but _NOR stands for NOV. A code with O
#              and 0 swapped is used only when the code itself has no clip.
#              BY_HAND assigns the clips whose name is wrong or cannot say (a
#              misspelled code, no code, two crops), with the evidence. All of
#              these are listed when the script runs. A recording without
#              exactly one clip keeps an empty video (the tool reports it
#              unreadable) and is listed with the unused clips of its folder and
#              trial whose code is at most two characters off. The frames
#              check finds a clip with other frames than its tracking rows,
#              but most clips have exactly 9000 or 18000 frames: only the
#              sync check (the mouse where the tracking puts it) tells two
#              animals' clips apart.
#   raw_dir    the test folder, where the clip check looks for the raw camera
#              recordings (GX*.MP4 in camera folders C1-C4; batch 1 plus maze
#              loose in the folder)
#   expected_s protocol length: 10 min social preference and open field,
#              5 min novel object and plus maze
#   arena      social preference 44 x 24 cm, novel object and open field
#              49 x 49 cm (area of the corners); plus maze 60 cm from tl to bl
#              (as epm_all.yaml) and the twelve corners of the maze as outline
#   objects    novel object objL objR (both trials: the habituation trial holds
#              two copies of the familiar object). Not the social preference
#              cups: their lids, the colour of the floor, hide most of each cup
#              and its centre, where the points sit, so the check sees only the
#              wire arc beside the lid and flagged 153 of 348 trials on the
#              light alone (2026-11-06 run). 28_cup_geometry.R checks the cup
#              points from the tracking.
#   object_set test and trial: the objects check compares batches within it
#   group      the batch
#
# Rows are ordered by batch, test and animal, so that the three social
# preference trials of an animal, cut from one camera recording, follow each
# other (the clip check then reads that recording once).
#
#   Rscript 32_video_qc_manifest.R
#   python -I tools/video_qc validation/exp9_boris/results/video_qc_manifest.csv validation/exp9_boris/results/video_qc --checks frames,objects
#   (each release, about 30 min; then sync, about 1.5 h, and clip, several
#   hours, with --checks sync or --checks clip into the same folder)
#
# Writes results/video_qc_manifest.csv (local paths; not kept in git) and
# prints the recordings without a clip and the clips no recording uses. The
# tool writes results/video_qc/: video_qc.csv (one row per recording, the
# flags of every check run there), which 31_light_checks.R reads, a table per
# check, and stills/ with --stills. Runtime under a minute (it lists the clip
# folders on the share).
# ============================================================================

suppressMessages({ library(stringr); library(data.table) })
SCRIPTS <- local({ a <- commandArgs(trailingOnly = FALSE); p <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(if (length(p)) dirname(p[[1]]) else ".", winslash = "/") })
PROJ <- normalizePath(file.path(SCRIPTS, ".."), winslash = "/")
RES <- file.path(PROJ, "results")
STAGE <- Sys.getenv("EXP9_STAGING", "C:/Users/topohl/Documents/exp9_sleap_staging")
VIDEOS <- Sys.getenv("EXP9_VIDEOS", "//mdc-berlin.net/fs/RD_Hoernberg/Tobi/E9_SIS")

# the staged recordings, listed as 31_light_checks.R lists them
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
recordings[, id := paste(assay, phase, batch, Code)]
stopifnot(!anyDuplicated(recordings$id))

# each test's folder per batch, and what the checks need to know about the test
FOLDERS <- data.table(batch = rep(paste0("B", 1:6), each = 4), assay = rep(c("SocP", "OFT", "NOR", "EPM"), 6),
  folder = c("20221117_SocP", "20221119_OFT", "20221120_NOR", "20221121_EPM", "20230218_SocP", "20230220_OFT", "20230221_NOR", "20230222_EPM",
             "SocP", "OFT", "NOR", "EPM", "SocP", "OFT", "NOR", "EPM",
             "20231210_SocP", "20231212_OFT", "20231213_NOR", "20231214_EPM", "SocP", "OFT", "NOR", "05012024_EPM"))
FOLDERS[, clip_dir := ifelse(batch == "B2" & assay == "NOR", "NOR_B2_cropped", "cropped")]
TEST <- data.table(assay = c("SocP", "NOR", "OFT", "EPM"), expected_s = c(600, 300, 600, 300),
  arena_cm = c("44x24", "49x49", "49x49", "60"), corners = c("tl tr br bl", "tl tr br bl", "tl tr br bl", "tl bl"),
  outline = c(NA, NA, NA, "tl tr ctr rt rb cbr br bl cbl lb lt ctl"), objects = c(NA, "objL objR", NA, NA))
TRIALS <- c("HAB", "S1", "S2", "NOV")

# Clips whose name is wrong or cannot say which recording they belong to (path from the test folder). Checked
# 2026-11-06: the clip has as many frames as the tracking has rows, and the tracked body centre lies on the mouse (video
# darker than its median background there) in 16 of 16 frames spread over the clip, 0-2 of 9-12 with the track shifted
# by 4000 frames.
BY_HAND <- data.table(
  id = c("OFT OFT B1 F9L3", "OFT OFT B1 V2P0", "EPM EPM B2 T4W6",
         "SocP HAB B2 Q6Q6", "SocP S1 B2 Q6Q6", "SocP S2 B2 Q6Q6", "SocP S2 B4 B7C9", "EPM EPM B1 M8P3"),
  clip = c("cropped/GX010385cropped.mp4", "cropped/old/384cropped.mp4", "cropped/E9_SIS_B2_T4W6_1.mp4",
           "cropped/E9_SIS_B2_O6O6_HAB.mp4", "cropped/E9_SIS_B2_O6O6_S1.mp4", "cropped/E9_SIS_B2_O6O6_S2.mp4",
           "cropped/E9_SIS_B4_B7C2_S2.mp4", "cropped/E9_B1_EPM_GX010035_M8P1.mp4"),
  why = c("no code in the name; cut from C2/GX010385; 18127 frames = rows; 16/16 on the mouse",
          "no code in the name; the tracked crop was moved to old/: 17980 frames = rows, 16/16 on the mouse; GX010384cropped.mp4 has 18102 frames, 3/16",
          "two crops; the re-crop _1 was tracked: 8900 frames = rows, 16/16 on the mouse; the first crop has 9000 frames (8 s of empty maze first), 3/16",
          "named O6O6 (no such animal); 18000 frames = rows; 16/16 on the mouse",
          "named O6O6 (no such animal); 18000 frames = rows; 16/16 on the mouse",
          "named O6O6 (no such animal); 18000 frames = rows; 16/16 on the mouse",
          "named B7C2 (no such animal); 17982 frames = rows; 16/16 on the mouse",
          "named M8P1 (no such animal); 9000 frames = rows; 16/16 on the mouse"))

# ---- the clips ----------------------------------------------------------------------------------------------------------
clips <- FOLDERS[, .(clip = list.files(file.path(VIDEOS, batch, folder, clip_dir), "[.](mp4|avi|mov)$", ignore.case = TRUE)),
                 by = .(batch, assay, folder, clip_dir)]
clips[, path := file.path(VIDEOS, batch, folder, clip_dir, clip)]
parts <- strsplit(toupper(sub("[.][^.]+$", "", clips$clip)), "[^A-Z0-9]+")
clips[, trial := vapply(parts, function(p) paste(intersect(p, TRIALS), collapse = " "), "")]
clips[, nor_for_nov := assay == "NOR" & trial == "" & vapply(parts, function(p) "NOR" %in% p, TRUE)]
clips[nor_for_nov == TRUE, trial := "NOV"]
clips[, recrop := vapply(parts, function(p) any(grepl("^[0-9]{1,2}$", p)), TRUE)]
clips[, other_batch := vapply(seq_along(parts), function(k) { b <- grep("^B[0-9]$", parts[[k]], value = TRUE); length(b) > 0 && !batch[k] %in% b }, TRUE)]
cat(nrow(clips), "clips in", nrow(FOLDERS), "folders\n")

# a code with O and 0 swapped, in every combination
swapped <- function(code) {
  ch <- strsplit(code, "")[[1]]; i <- which(ch %in% c("O", "0"))
  if (!length(i)) return(character())
  g <- expand.grid(rep(list(c("O", "0")), length(i)), stringsAsFactors = FALSE)
  setdiff(apply(g, 1, function(v) { x <- ch; x[i] <- v; paste(x, collapse = "") }), code)
}
has_code <- function(k, codes) vapply(k, function(j) any(codes %in% parts[[j]]), TRUE)

recordings[, want := ifelse(assay %in% c("SocP", "NOR"), phase, "")]
same_trial <- function(r) which(clips$assay == r$assay & clips$batch == r$batch & clips$trial == r$want)
found <- rbindlist(lapply(seq_len(nrow(recordings)), function(i) {
  r <- recordings[i]
  k <- same_trial(r)
  hit <- k[has_code(k, toupper(r$Code))]; how <- "name"
  if (!length(hit)) { hit <- k[has_code(k, swapped(toupper(r$Code)))]; how <- "O/0 swapped" }
  data.table(id = r$id, n = length(hit), how = if (length(hit)) how else "",
             clip_row = if (length(hit) == 1) hit else NA_integer_, candidates = paste(clips$clip[hit], collapse = " | "))
}))
recordings <- merge(recordings, found, by = "id", sort = FALSE)
recordings[, video := ifelse(is.na(clip_row), NA_character_, clips$path[clip_row])]

# by hand: the evidence is in BY_HAND
fold <- FOLDERS[, .(assay, batch, folder)]
bh <- merge(BY_HAND, recordings[, .(id, assay, batch, candidates)], by = "id")
bh <- merge(bh, fold, by = c("assay", "batch"))
bh[, path := file.path(VIDEOS, batch, folder, clip)]
stale <- c(setdiff(BY_HAND$id, bh$id), bh[!file.exists(path), sprintf("%s -> %s", id, clip)])
if (length(stale)) stop("BY_HAND names a recording that is not staged or a clip that is not there (renamed or moved? ",
                        "then update or drop the entry): ", paste(stale, collapse = "; "), call. = FALSE)
recordings[bh, on = "id", `:=`(video = i.path, how = "by hand", n = 1L)]

# one clip, one recording
twice <- recordings[!is.na(video), .(ids = paste(id, collapse = ", "), .N), by = video][N > 1]
if (nrow(twice)) recordings[video %in% twice$video, `:=`(video = NA_character_, how = "clip shared")]

# ---- report ---------------------------------------------------------------------------------------------------------------
cat(sprintf("%d of %d recordings have a clip\n", sum(!is.na(recordings$video)), nrow(recordings)))
if (nrow(twice)) { cat("\nclips that more than one recording matches (left empty):\n"); print(as.data.frame(twice)) }
# for a recording without a clip: unused clips of its folder and trial whose code is at most two characters off (a typo)
near <- function(r) {
  k <- same_trial(r); k <- k[!clips$path[k] %in% recordings$video]
  d <- vapply(k, function(j) { p <- parts[[j]][nchar(parts[[j]]) == nchar(r$Code)]; if (length(p)) min(adist(toupper(r$Code), p)) else Inf }, 0)
  if (!any(d <= 2)) return("")
  paste0(clips$clip[k][d <= 2], " (", d[d <= 2], " off)", collapse = " | ")
}
none <- recordings[is.na(video)]
if (nrow(none)) {
  none[, candidates := ifelse(how == "clip shared", "a clip it shares (above)", ifelse(n > 1, candidates, "none"))]
  none[, near_names := vapply(seq_len(.N), function(i) near(none[i]), "")]
  cat("\nrecordings without exactly one clip (video left empty; the tool reports them unreadable):\n")
  print(as.data.frame(none[, .(id, candidates, near_names)]))
}
cat("\nassigned by hand:\n"); print(as.data.frame(bh[, .(id, clip, matched_by_name = candidates, why)]))
note <- function(title, x) if (nrow(x)) { cat("\n", title, ":\n", sep = ""); print(as.data.frame(x)) }
used <- clips[, path] %in% recordings$video
note("matched with O and 0 swapped in the code", recordings[how == "O/0 swapped", .(id, video = basename(video))])
note("novel object clips named _NOR, read as _NOV", clips[used & nor_for_nov == TRUE, .(batch, clip)])
note("re-crops (a number in the name) in use", clips[used & recrop == TRUE, .(batch, assay, clip)])
note("clips whose name carries another batch", clips[other_batch == TRUE, .(batch, folder, clip)])
unused <- clips[!path %in% recordings$video]
cat("\nclips no recording uses:", nrow(unused), "\n")
invisible(unused[, cat(strwrap(paste0(batch, "/", folder, "/", clip_dir, " (", .N, "): ", paste(clip, collapse = " ")), 120, exdent = 4), sep = "\n"),
       by = .(batch, folder, clip_dir)])

# ---- manifest -------------------------------------------------------------------------------------------------------------
manifest <- merge(recordings, TEST, by = "assay", sort = FALSE)
manifest <- merge(manifest, fold, by = c("assay", "batch"), sort = FALSE)
manifest[, `:=`(tracking = file, raw_dir = file.path(VIDEOS, batch, folder), body = "bodycentre", group = batch,
                object_set = ifelse(is.na(objects), NA_character_, paste(assay, phase)))]
manifest <- manifest[, .(id, assay, phase, batch, Code, tracking, video, raw_dir, expected_s, arena_cm, corners, outline, objects,
                         object_set, body, group)]
manifest[, trial_order := match(phase, c("HAB", "S1", "S2", "NOV", "OFT", "EPM"))]
setorder(manifest, batch, assay, Code, trial_order)
manifest[, trial_order := NULL]
fwrite(manifest, file.path(RES, "video_qc_manifest.csv"), na = "")
cat("\nwrote", file.path(RES, "video_qc_manifest.csv"), "with", nrow(manifest), "recordings\n")
