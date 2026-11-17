# The SocP tracking QC report (validation/exp9_boris/scripts/24_qc_report.R) on
# a two-recording fixture: the script is copied into a temporary project with
# its own config, so it reads nothing outside the temporary folder.

qc_report_script <- file.path(repo_root, "validation", "exp9_boris", "scripts", "24_qc_report.R")

test_that("the QC report script parses", {
  expect_silent(parse(qc_report_script))
})

write_qc_report_fixture <- function(root) {
  csv <- function(d, ...) {
    path <- file.path(root, ...)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(d, path, row.names = FALSE)
  }
  codes <- c("B4Y7", "Q5M7")
  dir.create(file.path(root, "proj", "scripts"), recursive = TRUE)
  file.copy(qc_report_script, file.path(root, "proj", "scripts"))
  csv(data.frame(file = c("Q5M7_HAB", "B4Y7_HAB", "B4Y7_HAB"), target = c("socr", "socr", "all"),
                 start_s = c(328.3, 461.7, 0), end_s = c(329, 466, 0.07),
                 decision = c("not_on_object", "not_hidden", "leave_out"),
                 reason = c("rearing at the wire", "mouse in view", "borrowed frames"), reviewer = "test"),
      "proj", "config", "qc_decisions_socp.csv")
  csv(data.frame(file = "B4Y7_S1", kind = "end", first_frame = 17000L, last_frame = NA_integer_,
                 evidence = "lid lifted", reviewer = "test"), "proj", "config", "hand_events_socp.csv")
  dir.create(file.path(root, "A", "Exp9_manuscript", "config"), recursive = TRUE)
  writeLines(c("typography: {family: sans, axis_text_pt: 6, axis_title_pt: 7, legend_text_pt: 6, legend_title_pt: 7, annotation_pt: 5}",
               "group: {CON: '#888888', RES: '#1F3D52', SUS: '#F2CA4E'}",
               "line: {axis_pt: 0.5, reference_pt: 0.5}",
               "heatmap: {tile: 0.9}"),
             file.path(root, "A", "Exp9_manuscript", "config", "manuscript_palette.yml"))

  contact <- data.frame(file = paste0(codes, "_HAB"), Code = codes, contactNovel = c(20, 30),
                        contactFamiliar = c(10, 25), distance = c(3000, 3500))
  csv(contact, "ref", "assay_summaries", "SocP", "B4", "HAB", "combined_output.csv")
  dir.create(file.path(root, "ref", "data"))
  utils::write.table(data.frame(Code = codes, ID = c("OR1", "OR2"), Sex = c("Female", "Male"), Phenotype = c("SUS", "RES")),
                     file.path(root, "ref", "data", "exp9_sleap_all_animals.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  out <- file.path(root, "socp", "B4", "HAB")
  csv(transform(contact, contactNovel = c(21, 30), valid_fraction = 0.99, qcPass = TRUE, landmarksFlaggedS = c(1, 2),
                skeletonCollapseS = c(0, 0.5), noseReconstructedS = c(0.4, 0), bodycentreReconstructedS = c(0.1, 0),
                noseAnchorErrorP95 = c(0.8, 0.9), onCupLeft = 0, onCupRight = c(2, 0.5), contactRemovedLeft = 0,
                contactRemovedRight = c(1, 0.2)), "socp", "B4", "HAB", "combined_output.csv")
  csv(data.frame(Code = rep(codes, each = 2), landmark = c("nose", "bodycentre"), outside_s = 0, flagged_s = c(1, 0, 2, 0),
                 collapse_s = 0, reconstructed_s = c(0.4, 0.1, 0, 0), still_missing_s = 0),
      "socp", "B4", "HAB", "landmark_repairs.csv")
  csv(data.frame(Code = "B4Y7", check = "skeleton_split", target = "nose", start_s = 10, duration_s = 0.3,
                 reconstructed_share = 1, start_frame = 301L, end_frame = 309L), "socp", "B4", "HAB", "landmark_repair_episodes.csv")
  csv(data.frame(Code = codes, target = "nose", rank = 1L, anchor_a = "leftEar", anchor_b = "rightEar", frames = 900L,
                 error_median = 0.3, error_p95 = c(0.8, 0.9)), "socp", "B4", "HAB", "landmark_repair_models.csv")
  # two on-cup episodes at the right cup: one a not_on_object decision covers, one inside a not_hidden window
  csv(data.frame(Code = c("Q5M7", "B4Y7"), object = "socr", start_s = c(328.4, 462), duration_s = c(0.5, 2),
                 decision = "on_object", median_apparent_size = 1.2, start_frame = c(9853L, 13861L),
                 end_frame = c(9867L, 13920L)), "socp", "B4", "HAB", "on_cup_episodes.csv")
}

test_that("the QC report annotates on-cup episodes with on-cup decisions only, and lists the hands", {
  root <- temp_run_dir("qc_report")
  write_qc_report_fixture(root)
  previous <- Sys.getenv("E9_ANALYSIS_ROOT", unset = NA_character_)
  Sys.setenv(E9_ANALYSIS_ROOT = file.path(root, "A"))       # the child process inherits it
  on.exit(if (is.na(previous)) Sys.unsetenv("E9_ANALYSIS_ROOT") else Sys.setenv(E9_ANALYSIS_ROOT = previous), add = TRUE)
  out <- file.path(root, "report")
  output <- run_r_script(file.path(root, "proj", "scripts", "24_qc_report.R"),
                         c(file.path(root, "socp"), file.path(root, "ref"), out))
  expect_script_succeeded(output, "24_qc_report.R")

  cups <- utils::read.csv(file.path(out, "tables", "on_cup_episodes.csv"), stringsAsFactors = FALSE)
  expect_equal(cups$reviewed[cups$Code == "Q5M7"], "not_on_object: rearing at the wire")
  expect_true(cups$reviewed[cups$Code == "B4Y7"] %in% c("", NA))
  hands <- utils::read.csv(file.path(out, "tables", "hand_events.csv"), stringsAsFactors = FALSE)
  expect_equal(hands$file, "B4Y7_S1")
  expect_true(all(c("Reviewed decisions", "Hand events") %in% openxlsx::getSheetNames(file.path(out, "qc_report.xlsx"))))
})
