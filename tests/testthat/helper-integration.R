# Shared helpers for end-to-end assay runs.

#' Write a DLC/SLEAP-format CSV from a named list of n x 2 coordinate matrices.
write_tracking_fixture <- function(path, points, n) {
  header_names <- c("scorer")
  header_coords <- c("bodyparts")
  body <- data.frame(frame = 0:(n - 1))
  for (nm in names(points)) {
    header_names <- c(header_names, rep(nm, 3))
    header_coords <- c(header_coords, "x", "y", "likelihood")
    body[[paste0(nm, "_x")]] <- points[[nm]][, "x"]
    body[[paste0(nm, "_y")]] <- points[[nm]][, "y"]
    body[[paste0(nm, "_l")]] <- rep(1, n)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    paste(rep("scorer", length(header_names)), collapse = ","),
    paste(header_names, collapse = ","),
    paste(header_coords, collapse = ","),
    apply(body, 1, function(r) paste(r, collapse = ","))
  ), path)
  path
}

#' Run one of the batch assay scripts as a subprocess with a given config.
run_assay_script <- function(script_name, config_path) {
  script <- file.path(repo_root, "02_SLEAPanalzyer", script_name)
  rscript <- file.path(
    R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  # system2(env = ) is ignored on Windows, so set it here and let the child
  # process inherit it.
  previous <- Sys.getenv("SLEAP_ANALYZER_CONFIG", unset = NA_character_)
  Sys.setenv(SLEAP_ANALYZER_CONFIG = config_path)
  on.exit({
    if (is.na(previous)) Sys.unsetenv("SLEAP_ANALYZER_CONFIG") else
      Sys.setenv(SLEAP_ANALYZER_CONFIG = previous)
  }, add = TRUE)
  suppressWarnings(system2(
    rscript, c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE
  ))
}

#' Fail with the captured output when a script exits non-zero.
expect_script_succeeded <- function(output, label) {
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    testthat::fail(paste0(
      label, " exited with status ", status, ":\n", paste(output, collapse = "\n")
    ))
  }
  invisible(TRUE)
}

temp_run_dir <- function(prefix) {
  file.path(tempdir(), paste0(prefix, "-", as.integer(stats::runif(1, 1, 1e9))))
}

#' Square arena corner landmarks in pixels.
arena_corner_points <- function(n, side_px) {
  list(
    tl = cbind(x = rep(0, n), y = rep(side_px, n)),
    tr = cbind(x = rep(side_px, n), y = rep(side_px, n)),
    br = cbind(x = rep(side_px, n), y = rep(0, n)),
    bl = cbind(x = rep(0, n), y = rep(0, n))
  )
}
