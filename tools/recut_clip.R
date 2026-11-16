#!/usr/bin/env Rscript
# Re-cut a clip from its raw recording by frame numbers, with the crop of the clip it replaces, check it, and log it.
#
# A clip is a crop of a stretch of one raw camera file. When a clip turns out to be cut wrong (it stops before the end of
# the trial, starts during the handling, or holds another trial), the stretch is cut again from the raw file: frames
# <first> to <last> (0-based, as video tools count them), cropped to the same rectangle, re-encoded at the raw file's frame
# rate (H.264, constant quality) without sound. The new clip is then tracked like any other (SLEAP) and brought in with
# tools/retrack_recording.py.
#
# Checks, all before the log row is written: the frame range lies inside the raw file and the crop inside its picture; the
# new clip has exactly last - first + 1 frames and the crop's size; its first and last frames are the raw frames <first>
# and <last> (mean absolute grey difference of small greyscale copies, below --max-diff, and smaller than to the
# neighbouring raw frames). Nothing is overwritten: an existing output stops the run.
#
# usage:
#   Rscript tools/recut_clip.R --raw <raw.MP4> --first 120 --last 9119 --crop 446,0,1080,1080 --out <new_clip.mp4>
#                              [--log <cuts.csv>] [--reason "..."] [--ffmpeg <ffmpeg>] [--ffprobe <ffprobe>] [--crf 17]
#   --crop x,y,width,height in pixels of the raw picture (the replaced clip's crop; tools/video_qc records it as
#   crop_x_px, crop_y_px, crop_w_px, crop_h_px in video_qc_clip.csv). ffmpeg and ffprobe: the options, else the FFMPEG /
#   FFPROBE environment variables, else the PATH.
# Example (Exp9, 00662 / R8I2 plus maze, clip stopped 39 s before the end of the trial):
#   Rscript tools/recut_clip.R --raw //mdc-berlin.net/fs/RD_Hoernberg/Tobi/E9_SIS/B6/05012024_EPM/C3/GX011312.MP4 \
#     --first 120 --last 9119 --crop 446,0,1080,1080 \
#     --out //mdc-berlin.net/fs/RD_Hoernberg/Tobi/E9_SIS/B6/05012024_EPM/cropped/E9_B6_R8I2_EPM_1.mp4

recut_args <- function(args) {
  get <- function(name, default = NULL) { i <- match(paste0("--", name), args)
    if (is.na(i)) return(default); if (i == length(args)) stop("--", name, " needs a value", call. = FALSE); args[i + 1] }
  a <- list(raw = get("raw"), first = get("first"), last = get("last"), crop = get("crop"), out = get("out"),
            log = get("log"), reason = get("reason", ""), crf = get("crf", "17"), max_diff = get("max-diff", "6"),
            ffmpeg = get("ffmpeg", Sys.getenv("FFMPEG", unname(Sys.which("ffmpeg")))),
            ffprobe = get("ffprobe", Sys.getenv("FFPROBE", unname(Sys.which("ffprobe")))))
  for (k in c("raw", "first", "last", "crop", "out")) if (is.null(a[[k]])) stop("--", k, " is required", call. = FALSE)
  for (k in c("ffmpeg", "ffprobe")) if (!nzchar(a[[k]]) || !file.exists(a[[k]])) stop(k, " not found: pass --", k, " or set ", toupper(k), call. = FALSE)
  a$first <- as.integer(a$first); a$last <- as.integer(a$last); a$crf <- as.integer(a$crf); a$max_diff <- as.numeric(a$max_diff)
  a$crop <- as.integer(strsplit(a$crop, ",", fixed = TRUE)[[1]])
  stopifnot("--crop is x,y,width,height" = length(a$crop) == 4 && !anyNA(a$crop), "--first <= --last" = !is.na(a$first) && !is.na(a$last) && a$first <= a$last,
            "--first >= 0" = a$first >= 0, "even crop size (4:2:0 video)" = a$crop[3] %% 2 == 0 && a$crop[4] %% 2 == 0)
  a
}

# Width, height, frame rate (as a fraction) and frame count of the first video stream.
video_info <- function(ffprobe, path, count = FALSE) {
  out <- system2(ffprobe, c("-v", "error", "-select_streams", "v:0", if (count) "-count_frames",
                            "-show_entries", paste0("stream=width,height,r_frame_rate,", if (count) "nb_read_frames" else "nb_frames"),
                            "-of", "default=noprint_wrappers=1", shQuote(path)), stdout = TRUE, stderr = TRUE)
  kv <- setNames(sub("^[^=]*=", "", out), sub("=.*$", "", out))
  n <- if (count) kv[["nb_read_frames"]] else kv[["nb_frames"]]
  list(width = as.integer(kv[["width"]]), height = as.integer(kv[["height"]]), rate = kv[["r_frame_rate"]], frames = suppressWarnings(as.integer(n)))
}

# A 64 x 64 greyscale copy of frame n (0-based) of a video, optionally cropped first, as numbers 0-255.
grey_frame <- function(ffmpeg, path, n, crop = NULL) {
  vf <- paste0("select=eq(n\\,", n, ")", if (!is.null(crop)) sprintf(",crop=%d:%d:%d:%d", crop[3], crop[4], crop[1], crop[2]), ",scale=64:64,format=gray")
  tmp <- tempfile(fileext = ".gray")
  st <- system2(ffmpeg, c("-hide_banner", "-loglevel", "error", "-i", shQuote(path), "-vf", shQuote(vf), "-vsync", "0", "-frames:v", "1",
                          "-f", "rawvideo", "-y", shQuote(tmp)))
  if (st != 0 || !file.exists(tmp) || file.size(tmp) != 64 * 64) stop("could not read frame ", n, " of ", path, call. = FALSE)
  as.integer(readBin(tmp, "raw", 64 * 64))
}

recut_clip <- function(a) {
  if (!file.exists(a$raw)) stop("raw recording not found: ", a$raw, call. = FALSE)
  if (file.exists(a$out)) stop("output exists, not overwritten: ", a$out, call. = FALSE)
  if (!dir.exists(dirname(a$out))) stop("output folder does not exist: ", dirname(a$out), call. = FALSE)
  ri <- video_info(a$ffprobe, a$raw)
  if (is.na(ri$frames)) ri <- video_info(a$ffprobe, a$raw, count = TRUE)
  stopifnot("frame range inside the raw recording" = a$last < ri$frames,
            "crop inside the raw picture" = a$crop[1] + a$crop[3] <= ri$width && a$crop[2] + a$crop[4] <= ri$height)
  vf <- sprintf("select=between(n\\,%d\\,%d),setpts=N/FRAME_RATE/TB,crop=%d:%d:%d:%d", a$first, a$last, a$crop[3], a$crop[4], a$crop[1], a$crop[2])
  st <- system2(a$ffmpeg, c("-hide_banner", "-loglevel", "error", "-i", shQuote(a$raw), "-vf", shQuote(vf), "-r", ri$rate, "-vsync", "cfr",
                            "-an", "-c:v", "libx264", "-preset", "medium", "-crf", a$crf, "-pix_fmt", "yuv420p", shQuote(a$out)))
  if (st != 0 || !file.exists(a$out)) stop("ffmpeg failed (exit ", st, ")", call. = FALSE)
  # ---- checks
  ci <- video_info(a$ffprobe, a$out, count = TRUE)
  n_expected <- a$last - a$first + 1L
  diff_to <- function(k_out, k_raw) mean(abs(grey_frame(a$ffmpeg, a$out, k_out) - grey_frame(a$ffmpeg, a$raw, k_raw, a$crop)))
  d_first <- diff_to(0, a$first); d_last <- diff_to(n_expected - 1L, a$last)
  d_first_next <- if (a$first + 3L < ri$frames) diff_to(0, a$first + 3L) else NA_real_
  d_last_prev <- if (a$last >= 3L) diff_to(n_expected - 1L, a$last - 3L) else NA_real_
  checks <- c(frames = ci$frames == n_expected, size = ci$width == a$crop[3] && ci$height == a$crop[4], rate = ci$rate == ri$rate,
              first_frame = d_first < a$max_diff && (is.na(d_first_next) || d_first <= d_first_next),
              last_frame = d_last < a$max_diff && (is.na(d_last_prev) || d_last <= d_last_prev))
  res <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), raw = a$raw, first_frame = a$first, last_frame = a$last,
                    frames = ci$frames, rate = ci$rate, crop = paste(a$crop, collapse = ","), out = a$out,
                    sha256 = unname(tools::sha256sum(a$out)), first_frame_diff = round(d_first, 2), last_frame_diff = round(d_last, 2),
                    checks = paste(names(checks)[checks], collapse = " "), reason = a$reason, stringsAsFactors = FALSE)
  if (!all(checks)) { file.rename(a$out, paste0(a$out, ".failed")); print(res)
    stop("checks failed: ", paste(names(checks)[!checks], collapse = ", "), "; the clip was renamed to ", basename(a$out), ".failed", call. = FALSE) }
  if (!is.null(a$log)) utils::write.table(res, a$log, sep = ",", row.names = FALSE, col.names = !file.exists(a$log), append = file.exists(a$log), qmethod = "double")
  res
}

if (sys.nframe() == 0L) {
  res <- recut_clip(recut_args(commandArgs(trailingOnly = TRUE)))
  cat(sprintf("cut %s: raw frames %d-%d (%d frames, %s fps), crop %s; first/last frame difference %.2f / %.2f; sha256 %s\n",
              basename(res$out), res$first_frame, res$last_frame, res$frames, res$rate, res$crop, res$first_frame_diff, res$last_frame_diff, substr(res$sha256, 1, 12)))
}
