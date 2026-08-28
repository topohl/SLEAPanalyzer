# Deterministic core loader and thin wrappers around existing repository I/O.

behavior_core_files <- function() {
  c(
    "validation.R", "units.R", "tracking_data.R", "geometry.R",
    "interpolation.R", "dyadic.R", "homography.R",
    "calibration.R", "zones.R", "events.R", "metrics.R", "qc.R",
    "config.R", "assay_config.R", "provenance.R"
  )
}

source_behavior_core <- function(core_dir, envir = parent.frame()) {
  files <- file.path(core_dir, behavior_core_files())
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Missing behavioral core file(s): ", paste(basename(missing), collapse = ", "))
  }
  invisible(lapply(files, sys.source, envir = envir))
}

read_tracking_csv <- function(file, fps, reader = NULL) {
  validate_fps(fps)
  if (is.null(reader)) {
    if (!exists("ReadDLCDataFromCSV", mode = "function", inherits = TRUE)) {
      stop("ReadDLCDataFromCSV must be sourced before read_tracking_csv()")
    }
    reader <- get("ReadDLCDataFromCSV", mode = "function", inherits = TRUE)
  }
  tracking <- reader(file = file, fps = fps)
  validate_tracking_data(tracking)
  tracking
}

empty_metadata_table <- function(columns) {
  as.data.frame(
    stats::setNames(replicate(length(columns), character(), simplify = FALSE), columns),
    stringsAsFactors = FALSE
  )
}

read_metadata_table <- function(path, required_columns, sep = "",
                                if_missing = c("error", "empty")) {
  if_missing <- match.arg(if_missing)
  if (!file.exists(path)) {
    if (if_missing == "error") stop("Metadata file not found: ", path)
    warning("Metadata file not found: ", path)
    return(empty_metadata_table(required_columns))
  }
  out <- utils::read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  missing <- setdiff(required_columns, names(out))
  if (length(missing) > 0) {
    stop("Metadata file is missing column(s): ", paste(missing, collapse = ", "))
  }
  out
}

metadata_lookup <- function(data, code, value_column, code_column = "Code") {
  required <- c(code_column, value_column)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Metadata table is missing column(s): ", paste(missing, collapse = ", "))
  }
  values <- data[data[[code_column]] == code, value_column]
  if (length(values) == 0) return(NA_character_)
  if (length(values) > 1) warning("Multiple metadata rows found for code ", code, "; using the first")
  as.character(values[1])
}

#' Row-bind report rows that may not share the same columns.
#'
#' Replaces an undeclared data.table::rbindlist() dependency. Missing columns
#' are filled with NA rather than dropped, so a file that produced fewer
#' report entries than the others still contributes a row and the gap is
#' visible.
#'
#' @param rows a list of one-row lists or data frames
#' @return a data frame with the union of all columns
bind_report_rows <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.frame())
  frames <- lapply(rows, function(row) {
    if (is.data.frame(row)) return(row)
    as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  })
  columns <- unique(unlist(lapply(frames, names), use.names = FALSE))
  filled <- lapply(frames, function(frame) {
    absent <- setdiff(columns, names(frame))
    for (column in absent) frame[[column]] <- NA
    frame[, columns, drop = FALSE]
  })
  do.call(rbind, c(filled, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
}
