# Generic validation helpers for the dependency-light behavioral core.

validate_scalar_number <- function(value, name, positive = FALSE, allow_zero = FALSE) {
  if (length(value) != 1 || !is.numeric(value) || !is.finite(value)) {
    stop(name, " must be one finite numeric value")
  }
  if (positive && (value < 0 || (!allow_zero && value == 0))) {
    qualifier <- if (allow_zero) "non-negative" else "positive"
    stop(name, " must be ", qualifier)
  }
  invisible(value)
}

validate_fps <- function(fps) {
  validate_scalar_number(fps, "fps", positive = TRUE)
  invisible(fps)
}

validate_nonnegative_integer <- function(value, name) {
  validate_scalar_number(value, name, positive = TRUE, allow_zero = TRUE)
  if (value != floor(value)) stop(name, " must be an integer")
  invisible(as.integer(value))
}

validate_positive_integer <- function(value, name) {
  validate_scalar_number(value, name, positive = TRUE)
  if (value != floor(value)) stop(name, " must be an integer")
  invisible(as.integer(value))
}

validate_equal_lengths <- function(..., names = NULL) {
  values <- list(...)
  lengths <- vapply(values, length, integer(1))
  if (length(unique(lengths)) > 1) {
    if (is.null(names)) names <- paste0("value", seq_along(values))
    stop("Lengths differ: ", paste(paste0(names, "=", lengths), collapse = ", "))
  }
  invisible(lengths[1])
}

validate_numeric_vector <- function(value, name, finite = FALSE) {
  if (!is.numeric(value)) stop(name, " must be numeric")
  if (finite && any(!is.finite(value))) stop(name, " must contain only finite values")
  invisible(value)
}
