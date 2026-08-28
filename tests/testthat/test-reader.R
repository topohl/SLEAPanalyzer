# ReadDLCDataFromCSV input validation.

write_dlc_csv <- function(frames, points = list(bodycentre = NULL), path = NULL) {
  n <- length(frames)
  if (is.null(path)) path <- tempfile(fileext = ".csv")
  names_row <- c("scorer")
  coords_row <- c("bodyparts")
  body <- data.frame(frame = frames)
  for (nm in names(points)) {
    xy <- points[[nm]]
    if (is.null(xy)) xy <- cbind(x = seq_len(n), y = seq_len(n))
    names_row <- c(names_row, rep(nm, 3))
    coords_row <- c(coords_row, "x", "y", "likelihood")
    body[[paste0(nm, "_x")]] <- xy[, 1]
    body[[paste0(nm, "_y")]] <- xy[, 2]
    body[[paste0(nm, "_l")]] <- rep(1, n)
  }
  lines <- c(
    paste(rep("scorer", length(names_row)), collapse = ","),
    paste(names_row, collapse = ","),
    paste(coords_row, collapse = ","),
    apply(body, 1, function(r) paste(r, collapse = ","))
  )
  writeLines(lines, path)
  path
}

test_that("contiguous frame numbering is accepted", {
  path <- write_dlc_csv(0:9)
  tracking <- ReadDLCDataFromCSV(path, fps = 30)
  expect_equal(length(tracking$frames), 10L)
  expect_equal(tracking$fps, 30)
  expect_equal(tracking$distance.units, "pixel")
})

test_that("skipped frames are rejected rather than treated as adjacent", {
  # Every displacement is computed between adjacent rows, so a gap in the
  # frame numbering would silently inflate speed and distance.
  path <- write_dlc_csv(c(0:4, 10:14))
  expect_error(ReadDLCDataFromCSV(path, fps = 30), "not contiguous")
  expect_error(ReadDLCDataFromCSV(path, fps = 30), "5 frame\\(s\\) skipped")
})

test_that("non-monotonic frame numbering is rejected", {
  path <- write_dlc_csv(c(0, 1, 3, 2, 4))
  expect_error(ReadDLCDataFromCSV(path, fps = 30), "monotonically")
})

test_that("duplicate frames are rejected", {
  path <- write_dlc_csv(c(0, 1, 1, 2))
  expect_error(ReadDLCDataFromCSV(path, fps = 30), "must be unique")
})

test_that("a constant stride greater than one warns about the sampled rate", {
  path <- write_dlc_csv(seq(0, 18, by = 2))
  expect_warning(ReadDLCDataFromCSV(path, fps = 15), "increments by 2")
})

test_that("an invalid fps is rejected", {
  path <- write_dlc_csv(0:9)
  expect_error(ReadDLCDataFromCSV(path, fps = 0), "positive finite")
  expect_error(ReadDLCDataFromCSV(path, fps = -1), "positive finite")
})

test_that("a missing file is reported clearly", {
  expect_error(
    ReadDLCDataFromCSV(file.path(tempdir(), "definitely-absent.csv"), fps = 30),
    "does not exist"
  )
})
