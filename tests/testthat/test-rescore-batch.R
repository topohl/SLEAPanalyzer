# Tests of tools/rescore_batch.R: the configuration copy, the backup, the run
# manifest check, the plan, and whole rescorings from a throwaway git checkout
# with stub assay scripts and one tiny real NOR run. No real data is touched.

source(file.path(repo_root, "tools", "compare_outputs.R"), local = TRUE)
source(file.path(repo_root, "tools", "rescore_batch.R"), local = TRUE)

run_git <- function(repo, ...) {
  out <- suppressWarnings(system2("git", c("-C", shQuote(repo), ...), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) stop("git ", paste(c(...), collapse = " "), ": ", paste(out, collapse = "\n"))
  out
}

write_text <- function(path, lines) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  path
}

commit_all <- function(repo, message = "fixture") {
  run_git(repo, "add", "-A")
  run_git(repo, "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-q", "-m", shQuote(message))
}

# A git checkout in root/checkout holding the given files (path -> lines), all committed.
make_checkout <- function(root, files, copy_code = FALSE) {
  repo <- file.path(root, "checkout")
  dir.create(repo, recursive = TRUE, showWarnings = FALSE)
  if (copy_code) file.copy(file.path(repo_root, "02_SLEAPanalzyer"), repo, recursive = TRUE)
  for (p in names(files)) write_text(file.path(repo, p), files[[p]])
  run_git(repo, "init", "-q")
  run_git(repo, "config", "core.autocrlf", "false")
  commit_all(repo)
  normalizePath(repo, winslash = "/")
}

new_root <- function(prefix) {
  root <- temp_run_dir(prefix)
  dir.create(root, recursive = TRUE)
  normalizePath(root, winslash = "/")
}

DECISIONS_HEADER <- "file,target,start_s,end_s,decision,reason,reviewer"
DAY <- as.Date("2026-11-16")

# Stub NOR scoring: distance is the number in each tracking file plus its leave-out rows.
stub_nor_script <- function(commit = NULL) c(
  "cfg <- yaml::read_yaml(Sys.getenv('SLEAP_ANALYZER_CONFIG'))",
  "decisions <- utils::read.csv(cfg$qc_decisions_file, stringsAsFactors = FALSE)",
  "files <- list.files(file.path(cfg$input_dir, cfg$batches), pattern = '[.]csv$', full.names = TRUE)",
  "codes <- tools::file_path_sans_ext(basename(files))",
  "distance <- vapply(files, function(f) as.numeric(readLines(f)[1]), 0) +",
  "  vapply(codes, function(code) sum(decisions$file == code), 0)",
  "out <- file.path(cfg$output_dir, cfg$batches)",
  "dir.create(out, recursive = TRUE, showWarnings = FALSE)",
  "utils::write.csv(data.frame(Code = codes, distance = unname(distance)), file.path(out, 'combined_output.csv'), row.names = FALSE)",
  "commit <- system2('git', c('rev-parse', 'HEAD'), stdout = TRUE)",
  if (!is.null(commit)) sprintf("commit <- '%s'", commit),
  "dirty <- length(system2('git', c('status', '--porcelain'), stdout = TRUE)) > 0",
  "writeLines(yaml::as.yaml(list(code = list(commit = commit, dirty_working_tree = dirty),",
  "                              extra = list(assay = 'NOR', batch = cfg$batches))), file.path(out, 'run_manifest.yaml'))"
)

# A stub NOR set-up: staged tracking of two animals, a live output of batch 1 and a checkout.
make_stub_nor <- function(prefix, script = stub_nor_script()) {
  root <- new_root(prefix)
  write_text(file.path(root, "stage", "NOR", "B1", "A001.csv"), "10")
  write_text(file.path(root, "stage", "NOR", "B1", "A002.csv"), "20")
  live <- file.path(root, "project", "sleap_output_all", "NOR", "B1")
  write_text(file.path(live, "combined_output.csv"), c("\"Code\",\"distance\"", "\"A001\",10", "\"A002\",20"))
  write_text(file.path(live, "plots", "A001.svg"), "<svg/>")
  write_text(file.path(live, "run_manifest.yaml"), "code: {commit: earlier}")
  repo <- make_checkout(root, list(
    "02_SLEAPanalzyer/DLCA_NOR v1.2.1.R" = script,
    "validation/exp9_boris/config/nor_all.yaml" = c(
      "# Stub NOR configuration.",
      paste0("input_dir: ", root, "/stage/NOR"),
      paste0("output_dir: ", root, "/project/sleap_output_all/NOR   # the live output"),
      "batches: [B1, B2]",
      "qc_decisions_file: qc_decisions_nor.csv  # resolved against this folder",
      "write_manifest: true"),
    "validation/exp9_boris/config/qc_decisions_nor.csv" = DECISIONS_HEADER
  ))
  list(root = root, repo = repo, live = live, decisions = file.path(repo, "validation/exp9_boris/config/qc_decisions_nor.csv"),
       backup_root = file.path(root, "project", "_superseded_2026-11-16_t1"))
}

add_decision <- function(f, code, end_s = 1) {
  cat(sprintf("%s,all,0,%s,leave_out,fixture,fixture\n", code, end_s), file = f$decisions, append = TRUE)
  commit_all(f$repo, paste("Leave out the start of", code))
}

superseded <- function(dir) list.files(dir, pattern = "^_superseded_")

test_that("the configuration copy sets batches and phases and makes relative paths absolute", {
  dir <- normalizePath(tempdir(), winslash = "/")
  lines <- c(
    "# A comment stays.",
    "input_dir: data/in   # relative",
    "output_dir: 'out dir/x'",
    "metadata_dir:",
    "animal_id_code_file: ~",
    "zone_file: null",
    "qc_decisions_file: \"dec isions.csv\" # quoted",
    "other: rel/not/a/path/key",
    "batches:",
    "  - B1",
    "  - B2",
    "phases: [S1, S2]  # both",
    "fps: 30"
  )
  set <- list(batches = "B2", phases = "S1")
  copy <- rewrite_config_lines(lines, set, rescore_path_keys("SocP"), dir)
  expect_equal(copy$lines, c(
    "# A comment stays.",
    paste0("input_dir: ", dir, "/data/in   # relative"),
    paste0("output_dir: '", dir, "/out dir/x'"),
    "metadata_dir:",
    "animal_id_code_file: ~",
    "zone_file: null",
    paste0("qc_decisions_file: \"", dir, "/dec isions.csv\" # quoted"),
    "other: rel/not/a/path/key",
    "batches: [B2]",
    "phases: [S1]",
    "fps: 30"
  ))
  expect_equal(copy$changes$key, c("batches", "phases", "input_dir", "output_dir", "qc_decisions_file"))
  expect_equal(check_config_copy(lines, copy$lines, set, rescore_path_keys("SocP"), dir), character())
  parsed <- yaml::yaml.load(paste(copy$lines, collapse = "\n"))
  expect_equal(parsed$qc_decisions_file, paste0(dir, "/dec isions.csv"))
  expect_null(parsed$animal_id_code_file)

  # Absolute paths stay, a missing batches key is added, and the OFT keys are its own.
  kept <- c("input_dir: C:/abs/x", "output_dir: \\\\server\\share\\x", "metadata_dir: ~/x", "behavior_root: ../stage/OFT")
  copy <- rewrite_config_lines(kept, list(batches = "B3"), rescore_path_keys("NOR"), dir)
  expect_equal(copy$lines, c(kept, "batches: [B3]"))
  copy <- rewrite_config_lines(kept, list(batches = "B3"), rescore_path_keys("OFT"), dir)
  expect_equal(copy$lines[4], paste0("behavior_root: ", normalizePath(file.path(dir, "../stage/OFT"), winslash = "/", mustWork = FALSE)))

  # The check notices a copy that differs anywhere else.
  bad <- sub("fps: 30", "fps: 25", rewrite_config_lines(lines, set, rescore_path_keys("SocP"), dir)$lines)
  expect_equal(check_config_copy(lines, bad, set, rescore_path_keys("SocP"), dir), "fps changed")
  expect_error(rewrite_config_lines(c("batches: [B1]", "batches: [B2]"), list(batches = "B1"), character(), dir), "more than once")
})

test_that("backups copy every file, record md5 sums and never overwrite", {
  root <- new_root("rescore-backup")
  src <- file.path(root, "out", "B1")
  write_text(file.path(src, "combined_output.csv"), c("Code,x", "A,1"))
  write_text(file.path(src, "plots", "p.svg"), "<svg/>")
  dst <- file.path(root, "_superseded_x", "NOR_B1_before")
  manifest <- file.path(root, "_superseded_x", "NOR_B1_before_md5.csv")
  dir.create(dirname(dst))
  backup_folder(src, dst, manifest)
  rows <- utils::read.csv(manifest, colClasses = "character")
  expect_equal(rows$file, c("combined_output.csv", "plots/p.svg"))
  expect_equal(rows$md5, unname(tools::md5sum(file.path(src, rows$file))))
  expect_equal(unname(tools::md5sum(file.path(dst, rows$file))), rows$md5)
  expect_error(backup_folder(src, dst, manifest), "will not overwrite")
  expect_equal(free_run_id(dirname(dst), "NOR_B1"), "NOR_B1__2")
  expect_equal(free_run_id(dirname(dst), "NOR_B2"), "NOR_B2")

  # No output yet: an empty backup and an empty manifest.
  backup_folder(file.path(root, "missing"), file.path(root, "_superseded_x", "NOR_B2_before"), file.path(root, "_superseded_x", "NOR_B2_before_md5.csv"))
  expect_true(dir.exists(file.path(root, "_superseded_x", "NOR_B2_before")))
  expect_equal(nrow(utils::read.csv(file.path(root, "_superseded_x", "NOR_B2_before_md5.csv"))), 0L)
})

test_that("the run manifest check wants this run, the checkout's HEAD, a clean tree, the batch and the phase", {
  root <- new_root("rescore-manifest")
  path <- write_text(file.path(root, "run_manifest.yaml"), yaml::as.yaml(list(
    code = list(commit = "abc123", dirty_working_tree = FALSE), extra = list(batch = "B1", phase = "S1"))))
  ok <- check_run_manifest(path, "abc123", "B1", "S1", since = Sys.time() - 60)
  expect_true(all(ok$ok))
  expect_equal(nrow(ok), 5L)
  expect_equal(check_run_manifest(path, "abc124", "B1")$ok, c(FALSE, TRUE, TRUE))
  expect_false(all(check_run_manifest(path, "abc123", "B1", "S2")$ok))
  expect_false(check_run_manifest(path, "abc123", "B1", since = Sys.time() + 60)$ok[1])
  write_text(path, yaml::as.yaml(list(code = list(commit = "abc123", dirty_working_tree = TRUE), extra = list(batch = "B1"))))
  expect_equal(check_run_manifest(path, "abc123", "B1")$ok, c(TRUE, FALSE, TRUE))
  missing <- check_run_manifest(file.path(root, "none.yaml"), "abc123", "B1")
  expect_equal(missing$ok, FALSE)
})

test_that("the plan names every path, and a plan only writes nothing", {
  f <- make_stub_nor("rescore-plan")
  plan <- rescore_plan("NOR", "B1", repo = f$repo, tag = "t1", date = DAY)
  expect_equal(plan$problems, character())
  expect_equal(c(plan$phase, plan$prefix, plan$id), c("NOV", "NOR/B1", "NOR_B1"))
  expect_equal(plan$head, run_git(f$repo, "rev-parse", "HEAD"))
  expect_equal(plan$out_dir, f$live)
  expect_equal(plan$backup_dir, file.path(f$backup_root, "NOR_B1_before"))
  expect_equal(plan$config_copy, file.path(f$backup_root, "NOR_B1.yaml"))
  expect_true("batches: [B1]" %in% plan$copy_lines)
  expect_true(paste0("qc_decisions_file: ", f$decisions, "  # resolved against this folder") %in% plan$copy_lines)
  printed <- capture.output(print_rescore_plan(plan))
  expect_true(any(grepl("qc_decisions_file: qc_decisions_nor.csv -> ", printed, fixed = TRUE)))

  capture.output(status <- rescore_batch_main(c("--assay", "NOR", "--batch", "B1", "--repo", f$repo, "--tag", "t1")))
  expect_equal(status, 0L)
  expect_length(superseded(file.path(f$root, "project")), 0L)

  expect_true(any(grepl("B3 is not among the batches", rescore_plan("NOR", "B3", repo = f$repo, date = DAY)$problems)))
  expect_error(rescore_plan("EPM", "B1", phase = "S1", repo = f$repo), "no phases")
  expect_error(rescore_plan("SocP", "B1", repo = f$repo), "needs --phase")
  expect_error(rescore_plan("NOR", "B1", phase = "S1", repo = f$repo), "must be one of")
  expect_error(rescore_plan("NOR", "B1", repo = f$repo, tag = "a/b"), "--tag")
})

test_that("a checkout with uncommitted or untracked files is refused and nothing is written", {
  f <- make_stub_nor("rescore-dirty")
  cat("# edited\n", file = file.path(f$repo, "02_SLEAPanalzyer", "DLCA_NOR v1.2.1.R"), append = TRUE)
  plan <- rescore_plan("NOR", "B1", repo = f$repo, tag = "t1", date = DAY)
  expect_true(any(grepl("uncommitted changes to tracked files", plan$problems)))
  expect_error(run_rescore(plan), "will not run")
  capture.output(status <- rescore_batch_main(c("--assay", "NOR", "--batch", "B1", "--repo", f$repo, "--write")))
  expect_equal(status, 1L)
  expect_length(superseded(file.path(f$root, "project")), 0L)

  run_git(f$repo, "checkout", "--", ".")
  write_text(file.path(f$repo, "notes.txt"), "scratch")
  plan <- rescore_plan("NOR", "B1", repo = f$repo, tag = "t1", date = DAY)
  expect_equal(length(plan$problems), 1L)
  expect_match(plan$problems, "untracked files")
})

test_that("a stub NOR rescoring backs up, runs, checks the manifest and compares", {
  f <- make_stub_nor("rescore-run")
  add_decision(f, "A001")
  allow <- write_text(file.path(f$root, "allow.csv"), c("file,key,column,reason", "NOR/B1/*,A001,,left out"))
  before <- unname(tools::md5sum(file.path(f$live, "combined_output.csv")))

  plan <- rescore_plan("NOR", "B1", repo = f$repo, tag = "t1", allow = allow, date = DAY)
  printed <- capture.output(status <- run_rescore(plan))
  expect_equal(status, 0L, info = paste(printed, collapse = "\n"))
  expect_equal(unname(tools::md5sum(file.path(plan$backup_dir, "combined_output.csv"))), before)
  expect_setequal(utils::read.csv(plan$backup_md5)$file, c("combined_output.csv", "plots/A001.svg", "run_manifest.yaml"))
  copy <- yaml::read_yaml(plan$config_copy)
  expect_equal(copy$batches, "B1")
  expect_equal(copy$qc_decisions_file, f$decisions)
  diffs <- utils::read.csv(plan$compare_out, colClasses = "character")
  expect_equal(unlist(diffs[1, c("file", "change", "key", "column", "old", "new", "expected")], use.names = FALSE),
               c("NOR/B1/combined_output.csv", "cell", "A001", "distance", "10", "11", "TRUE"))
  expect_equal(nrow(diffs), 1L)
  log <- readLines(plan$log)
  expect_true(paste("commit:", plan$head) %in% log)
  expect_true("exit status: 0" %in% log)

  # A second rescoring with the same tag keeps the first backup; a change outside the allowlist fails it.
  add_decision(f, "A002")
  plan2 <- rescore_plan("NOR", "B1", repo = f$repo, tag = "t1", allow = allow, date = DAY)
  expect_equal(plan2$id, "NOR_B1__2")
  printed <- capture.output(status <- run_rescore(plan2))
  expect_equal(status, 1L)
  expect_true(any(grepl("UNEXPECTED NOR/B1/combined_output.csv cell A002 distance 20 -> 21", printed, fixed = TRUE)))
  expect_equal(unname(tools::md5sum(file.path(plan$backup_dir, "combined_output.csv"))), before)
})

test_that("a run manifest that does not name the checkout's HEAD fails the rescoring", {
  f <- make_stub_nor("rescore-wrong-commit", script = stub_nor_script(commit = "deadbeef"))
  plan <- rescore_plan("NOR", "B1", repo = f$repo, tag = "t1", date = DAY)
  printed <- capture.output(status <- run_rescore(plan))
  expect_equal(status, 1L)
  expect_true(any(grepl("FAILED    commit is the checkout's HEAD (deadbeef)", printed, fixed = TRUE)))
})

test_that("an OFT rescoring puts the all-batches folder back", {
  root <- new_root("rescore-oft")
  stage <- file.path(root, "stage", "OFT")
  write_text(file.path(stage, "B1", "X001.csv"), "6")
  write_text(file.path(stage, "B1", "output_v1.2.0", "combined_summary.csv"), c("Code,distance_cm", "X001,5"))
  global <- file.path(stage, "OFT_SLEAP_output_v1.2.0")
  write_text(file.path(global, "all_batches_summary.csv"), c("Batch,Code,distance_cm", "B1,X001,5", "B2,Y001,7"))
  global_md5 <- unname(tools::md5sum(file.path(global, "all_batches_summary.csv")))
  repo <- make_checkout(root, list(
    "02_SLEAPanalzyer/DLCA_OFT v1.2.0.R" = c(
      "cfg <- yaml::read_yaml(Sys.getenv('SLEAP_ANALYZER_CONFIG'))",
      "files <- list.files(file.path(cfg$behavior_root, cfg$batches, cfg$input_subdir), pattern = '[.]csv$', full.names = TRUE)",
      "d <- data.frame(Batch = cfg$batches, Code = tools::file_path_sans_ext(basename(files)),",
      "                distance_cm = vapply(files, function(f) as.numeric(readLines(f)[1]), 0))",
      "out <- file.path(cfg$behavior_root, cfg$batches, cfg$output_subdir)",
      "utils::write.csv(d[, -1], file.path(out, 'combined_summary.csv'), row.names = FALSE)",
      "g <- file.path(cfg$behavior_root, 'OFT_SLEAP_output_v1.2.0')",
      "utils::write.csv(d, file.path(g, 'all_batches_summary.csv'), row.names = FALSE)"),
    "validation/exp9_boris/config/oft_all.yaml" = c(
      "behavior_root: ../../../../stage/OFT   # relative to this folder",
      "input_subdir: \".\"",
      "output_subdir: output_v1.2.0",
      "batches: [B1, B2]")
  ))
  # OFT writes no manifest, so an untracked file does not stop it.
  write_text(file.path(repo, "notes.txt"), "scratch")

  plan <- rescore_plan("OFT", "B1", repo = repo, tag = "t1", date = DAY)
  expect_equal(plan$problems, character())
  expect_equal(plan$backup_root, file.path(root, "stage", "_superseded_2026-11-16_t1"))
  printed <- capture.output(status <- run_rescore(plan))
  expect_equal(status, 0L, info = paste(printed, collapse = "\n"))

  expect_equal(unname(tools::md5sum(file.path(global, "all_batches_summary.csv"))), global_md5)
  written <- utils::read.csv(file.path(plan$global_written, "all_batches_summary.csv"))
  expect_equal(written$Batch, "B1")
  expect_true(file.exists(file.path(plan$global_backup, "all_batches_summary.csv")))
  diffs <- utils::read.csv(plan$compare_out, colClasses = "character")
  expect_equal(unlist(diffs[1, c("file", "key", "column", "old", "new")], use.names = FALSE),
               c("OFT/B1/combined_summary.csv", "X001", "distance_cm", "5", "6"))
})

test_that("a tiny real NOR rescoring reads its relative decisions file and passes the manifest check", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("sp")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("stringr")

  root <- new_root("rescore-real-nor")
  n <- 300
  hold <- function(x, y) cbind(x = rep(x, n), y = rep(y, n))
  # 10 px per cm in a 49 x 49 cm arena; the nose rests on the left object for
  # the first `contact` frames, then in the middle of the arena.
  for (animal in list(c("A001", 60), c("A002", 30))) {
    contact <- as.integer(animal[2])
    nose <- cbind(x = c(rep(150, contact), rep(250, n - contact)), y = c(rep(215, contact), rep(100, n - contact)))
    write_tracking_fixture(file.path(root, "stage", "NOR", "B1", paste0(animal[1], ".csv")), c(
      arena_corner_points(n, 490),
      list(objL = hold(150, 245), objR = hold(350, 245), nose = nose,
           bodycentre = cbind(x = nose[, "x"], y = nose[, "y"] - 60))
    ), n)
  }
  write_text(file.path(root, "meta", "novelLoc.txt"), c("Code\tNovelLoc", "A001\tR", "A002\tR"))
  write_text(file.path(root, "meta", "animalIDCode.txt"), c("Code ID", "A001 mouse-1", "A002 mouse-2"))
  repo <- make_checkout(root, copy_code = TRUE, files = list(
    "validation/exp9_boris/config/nor_all.yaml" = c(
      paste0("input_dir: ", root, "/stage/NOR"),
      paste0("output_dir: ", root, "/project/sleap_output_all/NOR"),
      paste0("metadata_dir: ", root, "/meta"),
      paste0("animal_id_code_file: ", root, "/meta/animalIDCode.txt"),
      "novel_location_file: novelLoc.txt",
      "batches: [B1]",
      "fps: 30",
      "arena_width_cm: 49",
      "arena_height_cm: 49",
      "arena_corner_names: [tl, tr, br, bl]",
      "max_interpolation_gap_s: 0.2",
      "likelihood_cutoff: ~",
      "max_plausible_speed_cm_s: 200",
      "qc_min_valid_fraction: 0.8",
      "qc_max_longest_gap_s: 5.0",
      "qc_max_interpolated_fraction: 0.2",
      "movement_cutoff_cm_s: 5",
      "integration_period_frames: 5",
      "contact_geometry: radial",
      "contact_distance_cm: 4",
      "body_exclusion_distance_cm: 1",
      "contact_angle_deg: [70, 290]",
      "proximity_range_cm: [4, 8]",
      "proximity_angle_deg: [90, 270]",
      "qc_decisions_file: qc_decisions_nor.csv",
      "min_bout_s: 0.0",
      "max_gap_s: 0.0",
      "report_rearing: false",
      "write_manifest: true"),
    "validation/exp9_boris/config/qc_decisions_nor.csv" = DECISIONS_HEADER
  ))

  # The first scoring of the batch: there is no output to back up yet.
  plan <- rescore_plan("NOR", "B1", repo = repo, tag = "first", date = DAY)
  printed <- capture.output(status <- run_rescore(plan))
  expect_equal(status, 0L, info = paste(c(printed, readLines(plan$log)), collapse = "\n"))
  expect_true(any(grepl("ok        commit is the checkout's HEAD", printed, fixed = TRUE)))

  # A committed leave-out of A001's first 30 frames (0 to 0.99 s) changes A001 alone.
  f <- list(repo = repo, decisions = file.path(repo, "validation/exp9_boris/config/qc_decisions_nor.csv"))
  add_decision(f, "A001", end_s = 0.99)
  allow <- write_text(file.path(root, "allow.csv"), c("file,key,column", "NOR/B1/*,A001,"))
  plan <- rescore_plan("NOR", "B1", repo = repo, tag = "second", allow = allow, date = DAY)
  printed <- capture.output(status <- run_rescore(plan))
  expect_equal(status, 0L, info = paste(c(printed, readLines(plan$log)), collapse = "\n"))
  diffs <- utils::read.csv(plan$compare_out, colClasses = "character")
  expect_true(nrow(diffs) > 0)
  expect_true(all(startsWith(diffs$key, "A001")))
  contact <- diffs[diffs$file == "NOR/B1/combined_output.csv" & diffs$column == "contactLeft", ]
  expect_equal(as.numeric(c(contact$old, contact$new)), c(2, 1), tolerance = 1e-6)
})
