test_that("recut_clip cuts a frame range with a crop, checks it and logs it", {
  ffmpeg <- Sys.getenv("FFMPEG", unname(Sys.which("ffmpeg"))); ffprobe <- Sys.getenv("FFPROBE", unname(Sys.which("ffprobe")))
  testthat::skip_if(!nzchar(ffmpeg) || !file.exists(ffmpeg) || !nzchar(ffprobe) || !file.exists(ffprobe), "ffmpeg / ffprobe not available")
  source(file.path(Sys.getenv("SLEAP_ANALYZER_REPO_ROOT"), "tools", "recut_clip.R"), local = TRUE)
  dir <- tempfile("recut"); dir.create(dir)
  raw <- file.path(dir, "raw.mp4")   # 120 frames of a moving test picture, 320 x 240 at 29.97 frames/s
  st <- system2(ffmpeg, c("-hide_banner", "-loglevel", "error", "-f", "lavfi", "-i", shQuote("testsrc=size=320x240:rate=30000/1001"),
                          "-frames:v", "120", "-c:v", "libx264", "-pix_fmt", "yuv420p", shQuote(raw)))
  testthat::skip_if(st != 0, "ffmpeg could not write the test video")
  a <- recut_args(c("--raw", raw, "--first", "10", "--last", "49", "--crop", "40,20,160,120", "--out", file.path(dir, "clip_1.mp4"),
                    "--log", file.path(dir, "cuts.csv"), "--reason", "test", "--ffmpeg", ffmpeg, "--ffprobe", ffprobe))
  res <- recut_clip(a)
  expect_equal(res$frames, 40L)
  expect_equal(res$rate, "30000/1001")
  expect_true(grepl("first_frame", res$checks) && grepl("last_frame", res$checks))
  log <- utils::read.csv(file.path(dir, "cuts.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(log), 1L)
  expect_equal(log$sha256, unname(tools::sha256sum(file.path(dir, "clip_1.mp4"))))
  # nothing is overwritten; a range past the end and a crop outside the picture stop the run
  expect_error(recut_clip(a), "output exists")
  expect_error(recut_clip(modifyList(a, list(out = file.path(dir, "b.mp4"), last = 130L))), "frame range")
  expect_error(recut_clip(modifyList(a, list(out = file.path(dir, "c.mp4"), crop = c(200L, 0L, 160L, 120L)))), "crop inside")
  expect_error(recut_args(c("--raw", raw, "--first", "10", "--last", "49", "--crop", "1,1,161,120", "--out", "x.mp4",
                            "--ffmpeg", ffmpeg, "--ffprobe", ffprobe)), "even crop")
})
