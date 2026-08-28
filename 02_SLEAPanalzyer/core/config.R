# Declarative run configuration.
#
# The goal is that another researcher can reproduce a run without editing
# source code. Every input path, output path, acquisition rate, arena
# dimension, landmark name and threshold lives in a YAML file that is read,
# validated and recorded in the run manifest.

#' Recursively merge a user configuration over a set of defaults.
#'
#' Lists are merged key by key so a config file only needs to state what it
#' changes. A NULL in the user config removes the key.
merge_config <- function(defaults, overrides) {
  if (is.null(overrides)) return(defaults)
  if (!is.list(defaults) || !is.list(overrides)) return(overrides)
  for (key in names(overrides)) {
    value <- overrides[[key]]
    if (is.list(value) && is.list(defaults[[key]]) &&
        !is.null(names(value)) && !is.null(names(defaults[[key]]))) {
      defaults[[key]] <- merge_config(defaults[[key]], value)
    } else {
      defaults[[key]] <- value
    }
  }
  defaults
}

#' Describe one configuration field.
#'
#' @param type one of "number", "integer", "string", "logical", "path",
#'   "number_vector", "string_vector"
#' @param required whether the field must be present
#' @param positive for numeric fields, require a value above zero
#' @param allow_zero with `positive`, permit zero as well, so the field is
#'   non-negative. Needed for thresholds where zero means "disabled", such as a
#'   minimum bout duration of zero seconds.
#' @param choices for string fields, the permitted values
config_field <- function(type, required = TRUE, positive = FALSE,
                         allow_zero = FALSE, choices = NULL, length = NULL) {
  list(type = type, required = required, positive = positive,
       allow_zero = allow_zero, choices = choices, length = length)
}

#' Validate a configuration against a schema.
#'
#' Every problem is collected and reported together, because fixing a config
#' file one error per run is needlessly slow.
#'
#' @param config the merged configuration
#' @param schema a named list of config_field() descriptions
#' @return the configuration, invisibly, or an error listing every problem
validate_config <- function(config, schema) {
  if (!is.list(config)) stop("config must be a list")
  problems <- character()

  for (name in names(schema)) {
    spec <- schema[[name]]
    value <- config[[name]]
    if (is.null(value)) {
      if (isTRUE(spec$required)) {
        problems <- c(problems, sprintf("'%s' is required but missing", name))
      }
      next
    }
    expected_length <- spec$length
    if (!is.null(expected_length) && length(value) != expected_length) {
      problems <- c(problems, sprintf(
        "'%s' must have length %d, got %d", name, expected_length, length(value)
      ))
      next
    }

    problem <- switch(
      spec$type,
      number = if (!is.numeric(value) || any(!is.finite(value))) {
        "must be a finite number"
      } else if (isTRUE(spec$positive) && any(if (isTRUE(spec$allow_zero)) value < 0 else value <= 0)) {
        if (isTRUE(spec$allow_zero)) "must be non-negative" else "must be positive"
      } else NULL,
      integer = if (!is.numeric(value) || any(!is.finite(value)) ||
                    any(value != floor(value))) {
        "must be a whole number"
      } else if (isTRUE(spec$positive) && any(if (isTRUE(spec$allow_zero)) value < 0 else value <= 0)) {
        if (isTRUE(spec$allow_zero)) "must be non-negative" else "must be positive"
      } else NULL,
      number_vector = if (!is.numeric(value) || any(!is.finite(value))) {
        "must be a numeric vector of finite values"
      } else if (isTRUE(spec$positive) && any(if (isTRUE(spec$allow_zero)) value < 0 else value <= 0)) {
        if (isTRUE(spec$allow_zero)) "must contain only non-negative values" else "must contain only positive values"
      } else NULL,
      string = if (!is.character(value) || length(value) != 1 || !nzchar(value)) {
        "must be one non-empty string"
      } else if (!is.null(spec$choices) && !value %in% spec$choices) {
        sprintf("must be one of: %s", paste(spec$choices, collapse = ", "))
      } else NULL,
      string_vector = if (!is.character(value) || any(!nzchar(value))) {
        "must be a character vector of non-empty strings"
      } else NULL,
      logical = if (!is.logical(value) || length(value) != 1 || is.na(value)) {
        "must be TRUE or FALSE"
      } else NULL,
      path = if (!is.character(value) || length(value) != 1 || !nzchar(value)) {
        "must be one non-empty path"
      } else NULL,
      stop("Unknown config field type: ", spec$type)
    )
    if (!is.null(problem)) {
      problems <- c(problems, sprintf("'%s' %s", name, problem))
    }
  }

  unknown <- setdiff(names(config), c(names(schema), ".config_path", ".config_dir"))
  if (length(unknown) > 0) {
    problems <- c(problems, sprintf(
      "unknown configuration key(s): %s", paste(unknown, collapse = ", ")
    ))
  }

  if (length(problems) > 0) {
    stop(
      "Invalid configuration:\n",
      paste0("  - ", problems, collapse = "\n")
    )
  }
  invisible(config)
}

#' Load, merge and validate a YAML run configuration.
#'
#' Relative paths in `path_fields` are resolved against the directory holding
#' the config file, so a config can be moved with its data.
#'
#' @param path path to a YAML configuration file
#' @param defaults a named list of defaults the file overrides
#' @param schema a named list of config_field() descriptions, or NULL
#' @param path_fields names of fields to resolve relative to the config file
load_analysis_config <- function(path, defaults = list(), schema = NULL,
                                 path_fields = character()) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Reading configuration files requires the yaml package.")
  }
  if (!file.exists(path)) stop("Configuration file not found: ", path)
  overrides <- yaml::yaml.load_file(path)
  if (is.null(overrides)) overrides <- list()
  if (!is.list(overrides)) stop("Configuration file must contain a mapping: ", path)

  config <- merge_config(defaults, overrides)
  config$.config_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  config$.config_dir <- dirname(config$.config_path)

  for (field in path_fields) {
    value <- config[[field]]
    if (is.null(value) || !is.character(value)) next
    absolute <- grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", value)
    config[[field]] <- ifelse(
      absolute, value, file.path(config$.config_dir, value)
    )
  }

  if (!is.null(schema)) validate_config(config, schema)
  config
}

#' Configuration values suitable for recording in a run manifest.
#'
#' Drops the internal bookkeeping keys but keeps the resolved config path.
config_for_manifest <- function(config) {
  recorded <- config[setdiff(names(config), ".config_dir")]
  recorded
}

#' Overlay a YAML configuration onto an existing config list.
#'
#' For scripts whose configuration has grown organically and is not yet
#' described by a schema, this removes the need to edit source code to run the
#' analysis elsewhere without forcing a full migration first. Keys absent from
#' the overlay keep their in-script value.
#'
#' Overriding a key with a value of a different type is refused, because it
#' almost always means a typo in the YAML rather than an intended change.
#'
#' @param config the configuration list defined in the script
#' @param path a YAML file, or NULL to read SLEAP_ANALYZER_CONFIG
#' @param path_fields keys to resolve relative to the config file directory
#' @param required whether a missing configuration is an error
apply_config_overlay <- function(config, path = NULL, path_fields = character(),
                                 required = FALSE) {
  if (is.null(path) || !nzchar(path)) path <- Sys.getenv("SLEAP_ANALYZER_CONFIG")
  if (!nzchar(path)) {
    if (required) stop("No configuration file supplied; set SLEAP_ANALYZER_CONFIG.")
    return(config)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Reading configuration files requires the yaml package.")
  }
  if (!file.exists(path)) stop("Configuration file not found: ", path)
  overrides <- yaml::yaml.load_file(path)
  if (is.null(overrides)) overrides <- list()
  if (!is.list(overrides)) stop("Configuration file must contain a mapping: ", path)

  unknown <- setdiff(names(overrides), names(config))
  if (length(unknown) > 0) {
    stop(
      "Configuration file sets key(s) this analysis does not use: ",
      paste(unknown, collapse = ", "),
      ". Check for a typo; the accepted keys are: ",
      paste(sort(names(config)), collapse = ", ")
    )
  }
  mismatched <- names(overrides)[vapply(names(overrides), function(key) {
    original <- config[[key]]
    replacement <- overrides[[key]]
    if (is.null(original) || is.null(replacement)) return(FALSE)
    !identical(is.numeric(original), is.numeric(replacement)) ||
      !identical(is.character(original), is.character(replacement)) ||
      !identical(is.logical(original), is.logical(replacement))
  }, logical(1))]
  if (length(mismatched) > 0) {
    stop(
      "Configuration file changes the type of key(s): ",
      paste(mismatched, collapse = ", ")
    )
  }

  merged <- merge_config(config, overrides)
  config_dir <- dirname(normalizePath(path, winslash = "/", mustWork = TRUE))
  for (field in path_fields) {
    value <- merged[[field]]
    if (is.null(value) || !is.character(value)) next
    absolute <- grepl("^(/|~|[A-Za-z]:[/\\]|\\\\)", value)
    merged[[field]] <- ifelse(absolute, value, file.path(config_dir, value))
  }
  merged$.config_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  merged
}
