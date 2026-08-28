# Experiment design description and validation.
#
# The statistical layer is deliberately separate from measurement extraction.
# Extraction produces behavioral measurements and QC; this layer decides how
# those measurements should be modelled, and it decides that from the design of
# the experiment, never from which packages happen to be installed.

#' Describe the design of an experiment.
#'
#' Naming these columns explicitly is what allows the model structure to be
#' derived rather than guessed. A pipeline that does not know which column
#' identifies an animal cannot know whether two rows are two animals or two
#' observations of one animal, and that distinction determines the entire
#' analysis.
#'
#' @param data one row per observation
#' @param response the measurement column to model
#' @param group the treatment or condition column of interest
#' @param animal_id column uniquely identifying an animal
#' @param sex,batch,cage,litter optional blocking or nuisance columns
#' @param within optional within-animal factor, for example assay phase or
#'   time bin. Its presence is what makes repeated measures possible.
#' @param covariates additional fixed-effect columns
experiment_design <- function(data, response, group, animal_id,
                              sex = NULL, batch = NULL, cage = NULL,
                              litter = NULL, within = NULL,
                              covariates = character()) {
  if (!is.data.frame(data)) stop("data must be a data frame")
  named <- c(
    response = response, group = group, animal_id = animal_id,
    sex = sex, batch = batch, cage = cage, litter = litter, within = within
  )
  named <- named[!vapply(named, is.null, logical(1))]
  all_columns <- c(unlist(named, use.names = FALSE), covariates)
  missing <- setdiff(all_columns, names(data))
  if (length(missing) > 0) {
    stop("data is missing column(s): ", paste(missing, collapse = ", "))
  }
  if (!is.numeric(data[[response]])) {
    stop("response column '", response, "' must be numeric")
  }

  structure(
    list(
      data = data,
      response = response,
      group = group,
      animal_id = animal_id,
      sex = sex,
      batch = batch,
      cage = cage,
      litter = litter,
      within = within,
      covariates = covariates
    ),
    class = "experiment_design"
  )
}

#' Observations per animal, the quantity that decides the model structure.
observations_per_animal <- function(design) {
  complete <- design$data[
    !is.na(design$data[[design$response]]) & !is.na(design$data[[design$group]]),
    , drop = FALSE
  ]
  if (nrow(complete) == 0) return(integer())
  table(complete[[design$animal_id]])
}

#' Check a design for the conditions that invalidate an analysis.
#'
#' Returns problems and warnings rather than throwing, so a caller can report
#' every issue with a dataset at once.
validate_experiment_design <- function(design) {
  stopifnot(inherits(design, "experiment_design"))
  data <- design$data
  problems <- character()
  warnings <- character()

  complete <- data[!is.na(data[[design$response]]) & !is.na(data[[design$group]]), , drop = FALSE]
  if (nrow(complete) < 3) {
    problems <- c(problems, "fewer than three complete observations")
  }

  group_values <- unique(complete[[design$group]])
  if (length(group_values) < 2) {
    problems <- c(problems, sprintf(
      "grouping column '%s' has fewer than two levels", design$group
    ))
  }

  per_animal <- observations_per_animal(design)
  if (length(per_animal) > 0 && any(per_animal > 1) && is.null(design$within)) {
    warnings <- c(warnings, paste(
      "some animals contribute multiple observations but no within-animal",
      "factor was declared. If these are time bins or repeated sessions, name",
      "them with `within`; otherwise they will be treated as exchangeable",
      "repeated measures."
    ))
  }

  # An animal must not appear in more than one treatment group unless the
  # design is genuinely crossover, which would need a within factor.
  if (length(per_animal) > 0) {
    groups_per_animal <- tapply(
      as.character(complete[[design$group]]), complete[[design$animal_id]],
      function(x) length(unique(x))
    )
    if (any(groups_per_animal > 1) && is.null(design$within)) {
      problems <- c(problems, sprintf(
        "some animals appear in more than one level of '%s' without a declared within-animal factor",
        design$group
      ))
    }
  }

  # Complete confounding between a blocking factor and the group makes their
  # effects unidentifiable.
  for (blocking in c("batch", "cage", "litter")) {
    column <- design[[blocking]]
    if (is.null(column)) next
    crossing <- table(complete[[column]], complete[[design$group]])
    if (nrow(crossing) > 1 && all(rowSums(crossing > 0) == 1)) {
      problems <- c(problems, sprintf(
        "'%s' is completely confounded with '%s': every %s contains only one group",
        column, design$group, blocking
      ))
    }
  }

  list(
    ok = length(problems) == 0,
    problems = problems,
    warnings = warnings,
    n_observations = nrow(complete),
    n_animals = length(per_animal),
    max_observations_per_animal = if (length(per_animal) == 0) 0L else max(per_animal),
    min_observations_per_animal = if (length(per_animal) == 0) 0L else min(per_animal)
  )
}

#' Decide the model structure from the design.
#'
#' The rule is simple and deterministic:
#'
#' * A random intercept for animal is included only when at least one animal
#'   contributes more than one observation. Fitting `(1|ID)` with a single
#'   observation per animal is not conservative, it is unidentifiable: lme4
#'   refuses it outright because the number of groups equals the number of
#'   observations.
#' * Blocking factors with more than one level enter as fixed effects.
#'
#' Crucially, this depends only on the data. It does not consult installed
#' packages. If a mixed model is required and the package is absent, that is an
#' error to be fixed, not a reason to silently fit a different model.
model_structure <- function(design) {
  stopifnot(inherits(design, "experiment_design"))
  per_animal <- observations_per_animal(design)
  repeated <- length(per_animal) > 0 && max(per_animal) > 1

  fixed <- c(design$group, design$covariates)
  for (blocking in c("sex", "batch", "cage", "litter")) {
    column <- design[[blocking]]
    if (is.null(column)) next
    if (length(unique(design$data[[column]])) > 1) fixed <- c(fixed, column)
  }
  if (!is.null(design$within) && repeated) fixed <- c(fixed, design$within)
  fixed <- unique(fixed)

  random <- if (repeated) paste0("(1|", design$animal_id, ")") else character()

  list(
    repeated_measures = repeated,
    type = if (repeated) "linear mixed model" else "linear model",
    fixed_effects = fixed,
    random_effects = random,
    formula = stats::as.formula(paste(
      design$response, "~", paste(c(fixed, random), collapse = " + ")
    )),
    rationale = if (repeated) {
      sprintf(
        "at least one animal contributes %d observations, so a random intercept for %s is identifiable",
        max(per_animal), design$animal_id
      )
    } else {
      sprintf(
        "every animal contributes exactly one observation, so a random intercept for %s is unidentifiable and is omitted",
        design$animal_id
      )
    }
  )
}
