# Run provenance.
#
# Every analysis run writes a manifest recording what code, what configuration
# and what inputs produced the outputs sitting next to it. Without this, a
# results table a year later cannot be traced to the version of the pipeline
# that made it, and the numbers cannot be reproduced or corrected.

OUTPUT_SCHEMA_VERSION <- "2.0.0"

#' Current repository commit, or NA outside a git checkout.
repository_commit <- function(repo_dir = ".") {
  result <- tryCatch(
    suppressWarnings(system2(
      "git", c("-C", shQuote(repo_dir), "rev-parse", "HEAD"),
      stdout = TRUE, stderr = FALSE
    )),
    error = function(e) NULL
  )
  if (is.null(result) || length(result) == 0 || !nzchar(result[1])) {
    return(NA_character_)
  }
  trimws(result[1])
}

#' Whether the working tree has uncommitted changes.
#'
#' A run made from a dirty tree cannot be reproduced from its commit alone, so
#' the manifest must say so.
repository_is_dirty <- function(repo_dir = ".") {
  result <- tryCatch(
    suppressWarnings(system2(
      "git", c("-C", shQuote(repo_dir), "status", "--porcelain"),
      stdout = TRUE, stderr = FALSE
    )),
    error = function(e) NULL
  )
  if (is.null(result)) return(NA)
  length(result) > 0
}

#' Installed versions of the packages a run depends on.
package_versions <- function(packages) {
  if (length(packages) == 0) return(list())
  versions <- lapply(packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(pkg))
  })
  stats::setNames(versions, packages)
}

#' Identify input files by path, size, modification time and MD5.
#'
#' The hash is what makes a rerun verifiable: a file can be moved or renamed
#' and still be recognised, and a silently edited file is detected.
input_file_records <- function(paths) {
  if (length(paths) == 0) return(list())
  lapply(paths, function(path) {
    if (!file.exists(path)) {
      return(list(path = path, exists = FALSE))
    }
    info <- file.info(path)
    list(
      path = normalizePath(path, winslash = "/", mustWork = FALSE),
      exists = TRUE,
      bytes = as.numeric(info$size),
      modified = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
      md5 = unname(tools::md5sum(path))
    )
  })
}

#' Build a run manifest.
#'
#' @param config the validated run configuration
#' @param inputs paths of the input files consumed by the run
#' @param packages names of packages whose versions should be recorded
#' @param repo_dir the repository to read the commit from
#' @param extra any additional named values to record
run_manifest <- function(config = list(), inputs = character(),
                         packages = character(), repo_dir = ".",
                         extra = list()) {
  list(
    schema = list(
      output_schema_version = OUTPUT_SCHEMA_VERSION,
      manifest_version = "1.0.0"
    ),
    run = list(
      timestamp_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
      user = unname(Sys.info()[["user"]]),
      hostname = unname(Sys.info()[["nodename"]])
    ),
    code = list(
      commit = repository_commit(repo_dir),
      dirty_working_tree = repository_is_dirty(repo_dir)
    ),
    environment = list(
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      platform = R.version$platform,
      packages = package_versions(packages)
    ),
    config = config_for_manifest(config),
    inputs = input_file_records(inputs),
    extra = extra
  )
}

#' Write a run manifest next to the outputs it describes.
#'
#' YAML is used so the manifest stays readable without tooling.
write_run_manifest <- function(manifest, path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Writing a run manifest requires the yaml package.")
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(yaml::as.yaml(manifest, indent.mapping.sequence = TRUE), path)
  invisible(path)
}

#' Warn loudly when a run cannot be reproduced from its recorded commit.
assert_reproducible_run <- function(manifest, strict = FALSE) {
  problems <- character()
  if (is.na(manifest$code$commit)) {
    problems <- c(problems, "the analysis is not running from a git checkout")
  }
  if (isTRUE(manifest$code$dirty_working_tree)) {
    problems <- c(problems, "the working tree has uncommitted changes")
  }
  if (length(problems) > 0) {
    message <- paste0(
      "This run cannot be reproduced from its recorded commit alone: ",
      paste(problems, collapse = "; "), "."
    )
    if (strict) stop(message) else warning(message)
  }
  invisible(length(problems) == 0)
}
