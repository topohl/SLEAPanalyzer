# Tests of tools/compare_outputs.R on small synthetic output trees.

source(file.path(repo_root, "tools", "compare_outputs.R"), local = TRUE)

# Write text byte for byte, so the tests control line endings.
write_bytes <- function(path, text) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeBin(charToRaw(text), path)
  path
}

status_of <- function(result, file) result$files$status[result$files$file == file]

test_that("identical files, equal cells and changed cells are told apart", {
  root <- temp_run_dir("compare-trees")
  old <- file.path(root, "old")
  new <- file.path(root, "new")
  write_bytes(file.path(old, "a.csv"), "Code,x\nA1,1\n")
  write_bytes(file.path(new, "a.csv"), "Code,x\nA1,1\n")
  write_bytes(file.path(old, "b.csv"), "Code,x\nA1,1.0\n")
  write_bytes(file.path(new, "b.csv"), "\"Code\",\"x\"\r\nA1,1\r\n")
  # Rows and columns in another order; one cell moves beyond the tolerance, one within it.
  write_bytes(file.path(old, "sub", "c.csv"), "Code,x,y\nA1,1,a\nB2,2,b\n")
  write_bytes(file.path(new, "sub", "c.csv"), "Code,y,x\nB2,b,2.5\nA1,a,1.0000000000001\n")
  write_bytes(file.path(old, "gone.csv"), "Code,x\nA1,1\n")
  write_bytes(file.path(new, "added.csv"), "Code,x\nA1,1\n")
  write_bytes(file.path(new, "notes.txt"), "not a table\n")

  result <- compare_output_trees(old, new)
  expect_setequal(result$files$file, c("a.csv", "b.csv", "sub/c.csv", "gone.csv", "added.csv"))
  expect_equal(status_of(result, "a.csv"), "identical")
  expect_equal(status_of(result, "b.csv"), "equal_cells")
  expect_equal(status_of(result, "sub/c.csv"), "changed")
  expect_equal(status_of(result, "gone.csv"), "only_old")
  expect_equal(status_of(result, "added.csv"), "only_new")

  cells <- result$diffs[result$diffs$change == "cell", ]
  expect_equal(nrow(cells), 1L)
  expect_equal(unlist(cells[1, c("file", "key", "column", "old", "new")], use.names = FALSE),
               c("sub/c.csv", "B2", "x", "2", "2.5"))
  expect_setequal(result$diffs$change, c("cell", "file_removed", "file_added"))
  expect_true(all(is.na(result$diffs$expected)))
  expect_equal(result$unexpected, 0L)

  prefixed <- compare_output_trees(old, new, prefix = "NOR/B1")
  expect_true(all(startsWith(prefixed$files$file, "NOR/B1/")))
})

test_that("added and removed columns and rows are reported, and keys fall back to file and the first column", {
  root <- temp_run_dir("compare-shape")
  write_bytes(file.path(root, "old", "qc.csv"), "file,v\nf1,1\nf2,2\n")
  write_bytes(file.path(root, "new", "qc.csv"), "file,v,w\nf1,1,9\nf3,3,9\n")
  write_bytes(file.path(root, "old", "named.csv"), "name,v,gone\nn1,1,x\n")
  write_bytes(file.path(root, "new", "named.csv"), "name,v\nn1,2\n")

  result <- compare_output_trees(file.path(root, "old"), file.path(root, "new"))
  qc <- result$diffs[result$diffs$file == "qc.csv", ]
  expect_equal(qc$change[qc$change != "cell"], c("column_added", "row_added", "row_removed"))
  expect_equal(qc$column[qc$change == "column_added"], "w")
  expect_equal(qc$key[qc$change == "row_added"], "f3")
  expect_equal(qc$key[qc$change == "row_removed"], "f2")
  row <- result$files[result$files$file == "qc.csv", ]
  expect_equal(row$key, "file")
  expect_equal(c(row$columns_added, row$rows_added, row$rows_removed, row$cells_changed), c(1L, 1L, 1L, 0L))

  named <- result$diffs[result$diffs$file == "named.csv", ]
  expect_equal(result$files$key[result$files$file == "named.csv"], "name")
  expect_equal(named$change, c("column_removed", "cell"))
  expect_equal(unlist(named[2, c("key", "column", "old", "new")], use.names = FALSE), c("n1", "v", "1", "2"))
})

test_that("repeated keys get occurrence numbers and multi-row tables their own keys", {
  root <- temp_run_dir("compare-keys")
  write_bytes(file.path(root, "old", "episodes.csv"), "file,start\nf1,10\nf1,20\n")
  write_bytes(file.path(root, "new", "episodes.csv"), "file,start\nf1,10\nf1,25\n")
  write_bytes(file.path(root, "old", "combined_bins.csv"), "file,Code,bin_s,d\nx,A1,60,1\nx,A1,120,2\n")
  write_bytes(file.path(root, "new", "combined_bins.csv"), "file,Code,bin_s,d\nx,A1,120,2\nx,A1,60,1.5\n")

  result <- compare_output_trees(file.path(root, "old"), file.path(root, "new"))
  ep <- result$diffs[result$diffs$file == "episodes.csv", ]
  expect_equal(unlist(ep[1, c("change", "key", "old", "new")], use.names = FALSE), c("cell", "f1 #2", "20", "25"))
  bins <- result$diffs[result$diffs$file == "combined_bins.csv", ]
  expect_equal(unlist(bins[1, c("key", "column", "old", "new")], use.names = FALSE), c("A1 | 60", "d", "1", "1.5"))
  expect_equal(nrow(bins), 1L)

  keyed <- compare_output_trees(file.path(root, "old"), file.path(root, "new"),
                                keys = list("episodes.csv" = c("file", "start")))
  ep <- keyed$diffs[keyed$diffs$file == "episodes.csv", ]
  expect_equal(ep$change, c("row_added", "row_removed"))
  expect_equal(ep$key, c("f1 | 25", "f1 | 20"))
})

test_that("numbers are equal within a tolerance relative to their size", {
  x <- c("1000000", "1", "NA", "Inf", "abc", "1", "", "TRUE")
  y <- c("1000000.0005", "1.0000000001", "NA", "Inf", "abd", "1.1", "", "TRUE")
  expect_equal(cells_equal(x, y, 1e-9), c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  expect_true(cells_equal("1", "1.1", 0.2))
  expect_false(cells_equal("Inf", "-Inf", 1))
  expect_false(cells_equal("NA", "0", 1))
})

test_that("the allowlist marks differences outside it as unexpected", {
  root <- temp_run_dir("compare-allow")
  old <- file.path(root, "old")
  new <- file.path(root, "new")
  write_bytes(file.path(old, "NOR", "B1", "combined_output.csv"), "file,Code,distance,contact\nA1L9_NOV,A1L9,10,1\nK4T2_NOV,K4T2,20,2\n")
  write_bytes(file.path(new, "NOR", "B1", "combined_output.csv"),
              "file,Code,distance,contact,leftOutS\nA1L9_NOV,A1L9,11,1,0.07\nK4T2_NOV,K4T2,20,3,0\n")
  write_bytes(file.path(old, "NOR", "B1", "tracking_qc.csv"), "file,frames\nA1L9_NOV,9000\nK4T2_NOV,9000\n")
  write_bytes(file.path(new, "NOR", "B1", "tracking_qc.csv"), "file,frames\nA1L9_NOV,8998\nK4T2_NOV,9000\n")
  allow <- write_bytes(file.path(root, "allow.csv"), paste0(
    "file,key,column,reason\n",
    "NOR/B1/*,A1L9,,the borrowed frames are left out\n",
    "combined_output.csv,,leftOut*,new bookkeeping column\n"))

  result <- compare_output_trees(old, new, allow = allow)
  bad <- result$diffs[!result$diffs$expected, ]
  expect_equal(result$unexpected, 1L)
  expect_equal(unlist(bad[1, c("file", "key", "column", "old", "new")], use.names = FALSE),
               c("NOR/B1/combined_output.csv", "K4T2", "contact", "2", "3"))
  expect_equal(result$files$unexpected, c(1L, 0L))

  # A row that names keys does not allow an added column, and one that names a column no added row.
  narrow <- data.frame(file = "*", key = "A1L9", column = "*", stringsAsFactors = FALSE)
  d <- data.frame(file = "NOR/B1/x.csv", change = c("column_added", "row_added", "cell", "cell"),
                  key = c("", "A1L9", "A1L9 | 60", "A1L91"), column = c("leftOutS", "", "d", "d"),
                  old = "", new = "", stringsAsFactors = FALSE)
  expect_equal(allowed_differences(d, narrow), c(FALSE, FALSE, TRUE, FALSE))
  wildcard <- data.frame(file = "x.csv", key = "A1*;K4T2", column = "", stringsAsFactors = FALSE)
  expect_equal(allowed_differences(d, wildcard), c(FALSE, TRUE, TRUE, TRUE))
})

test_that("release mode maps the release onto the assay output root and the staging tree", {
  root <- temp_run_dir("compare-release")
  release <- file.path(root, "releases", "2026-11-05_test")
  summaries <- file.path(release, "assay_summaries")
  output_root <- file.path(root, "sleap_output_all")
  stage_root <- file.path(root, "staging")
  write_bytes(file.path(summaries, "NOR", "B1", "combined_output.csv"), "Code,x\nA1,1\n")
  write_bytes(file.path(output_root, "NOR", "B1", "combined_output.csv"), "Code,x\nA1,2\n")
  write_bytes(file.path(summaries, "OFT", "B2", "combined_summary.csv"), "Code,x\nB1,1\n")
  write_bytes(file.path(stage_root, "OFT", "B2", "output_v1.2.0", "combined_summary.csv"), "Code,x\nB1,1\n")
  write_bytes(file.path(summaries, "SocP", "B1", "S1", "tracking_qc.csv"), "file,x\nf,1\n")
  write_bytes(file.path(summaries, "SocP", "B1", "S1", "run_manifest.yaml"), "code: {}\n")

  pairs <- release_pairs(release, output_root, stage_root)
  expect_equal(pairs$file, c("NOR/B1/combined_output.csv", "OFT/B2/combined_summary.csv", "SocP/B1/S1/tracking_qc.csv"))
  expect_equal(pairs$new[2], file.path(stage_root, "OFT", "B2", "output_v1.2.0", "combined_summary.csv"))
  expect_true(is.na(pairs$new[3]))

  result <- compare_release(summaries, output_root, stage_root)
  expect_equal(result$files$status, c("changed", "identical", "only_old"))

  write_bytes(file.path(summaries, "XYZ", "B1", "a.csv"), "Code\nA\n")
  expect_error(release_pairs(release, output_root, stage_root), "XYZ")
})

test_that("the command line writes the differences and sets the exit status", {
  root <- temp_run_dir("compare-cli")
  write_bytes(file.path(root, "old", "a.csv"), "Code,x\nA1,1\nB2,1\n")
  write_bytes(file.path(root, "new", "a.csv"), "Code,x\nA1,2\nB2,3\n")
  allow <- write_bytes(file.path(root, "allow.csv"), "file,key,column\na.csv,A1,x\n")
  out <- file.path(root, "diffs.csv")

  capture.output(status <- compare_outputs_main(c(file.path(root, "old"), file.path(root, "new"), "--out", out)))
  expect_equal(status, 0L)
  written <- utils::read.csv(out, colClasses = "character")
  expect_equal(names(written), c("file", "change", "key", "column", "old", "new", "expected"))
  expect_equal(written$key, c("A1", "B2"))

  capture.output(status <- compare_outputs_main(c(file.path(root, "old"), file.path(root, "new"), "--allow", allow)))
  expect_equal(status, 1L)
  capture.output(status <- compare_outputs_main(c("only-one-folder")))
  expect_equal(status, 2L)
  expect_equal(suppressMessages(compare_outputs_main(c("--nonsense"))), 2L)

  # Run as a script: the exit status reaches the shell.
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  output <- suppressWarnings(system2(rscript, c("--vanilla", shQuote(file.path(repo_root, "tools", "compare_outputs.R")),
                                                shQuote(file.path(root, "old")), shQuote(file.path(root, "new")),
                                                "--allow", shQuote(allow)), stdout = TRUE, stderr = TRUE))
  expect_equal(attr(output, "status"), 1L)
  expect_true(any(grepl("UNEXPECTED a.csv cell B2 x 1 -> 3", output, fixed = TRUE)))
})
