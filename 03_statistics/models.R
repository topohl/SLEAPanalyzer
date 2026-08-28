# Model fitting and multiplicity handling for the statistical layer.
#
# Two rules govern this file.
#
# 1. The model comes from the design, never from the package library. If a
#    design requires a mixed model and lmerTest is absent, that is an error.
#    Silently fitting a fixed-effects model instead changes the inference
#    without telling anyone.
#
# 2. Nothing here interprets a p-value. Results carry the adjusted value, the
#    family it was adjusted within, and the size of that family, so the reader
#    can judge the evidence.

#' Fit the model implied by a design.
#'
#' @param design an experiment_design
#' @return a list with the fit, the structure used and the design validation
fit_group_model <- function(design) {
  stopifnot(inherits(design, "experiment_design"))
  validation <- validate_experiment_design(design)
  if (!validation$ok) {
    stop(
      "The design does not support an analysis:\n",
      paste0("  - ", validation$problems, collapse = "\n")
    )
  }
  for (message_text in validation$warnings) warning(message_text)

  structure_spec <- model_structure(design)
  complete <- design$data[
    !is.na(design$data[[design$response]]) & !is.na(design$data[[design$group]]),
    , drop = FALSE
  ]

  if (structure_spec$repeated_measures) {
    if (!requireNamespace("lmerTest", quietly = TRUE)) {
      stop(
        "This design has repeated observations per animal and therefore ",
        "requires a mixed model, but the lmerTest package is not installed. ",
        "Install it rather than falling back to a model that ignores the ",
        "repeated-measures structure."
      )
    }
    fit <- lmerTest::lmer(structure_spec$formula, data = complete)
  } else {
    fit <- stats::lm(structure_spec$formula, data = complete)
  }

  list(
    fit = fit,
    structure = structure_spec,
    validation = validation,
    design = design,
    n_used = nrow(complete)
  )
}

#' Test the group term of a fitted model.
#'
#' @return a one-row data frame; `p_raw` is unadjusted
test_group_effect <- function(model) {
  design <- model$design
  table_out <- tryCatch(as.data.frame(stats::anova(model$fit)), error = function(e) NULL)
  if (is.null(table_out) || !design$group %in% rownames(table_out)) {
    return(data.frame(
      response = design$response,
      term = design$group,
      model = model$structure$type,
      statistic = NA_real_,
      df = NA_real_,
      p_raw = NA_real_,
      n_observations = model$n_used,
      n_animals = model$validation$n_animals,
      note = "the group term could not be extracted from the model",
      stringsAsFactors = FALSE
    ))
  }
  row <- table_out[design$group, , drop = FALSE]
  statistic_column <- intersect(c("F value", "F.value"), names(row))
  p_column <- intersect(c("Pr(>F)", "Pr..F."), names(row))
  df_column <- intersect(c("NumDF", "Df"), names(row))

  data.frame(
    response = design$response,
    term = design$group,
    model = model$structure$type,
    statistic = if (length(statistic_column)) row[[statistic_column[1]]] else NA_real_,
    df = if (length(df_column)) row[[df_column[1]]] else NA_real_,
    p_raw = if (length(p_column)) row[[p_column[1]]] else NA_real_,
    n_observations = model$n_used,
    n_animals = model$validation$n_animals,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Adjust p-values within an explicitly named family.
#'
#' The family is recorded on every row, along with its size, because an
#' adjusted p-value is meaningless without knowing what it was adjusted
#' across. Rows whose model failed are excluded from the family size rather
#' than inflating it, which is what `p.adjust()` would do if handed NAs.
#'
#' @param results a data frame with a p_raw column
#' @param method a method accepted by stats::p.adjust
#' @param family_label a short description of the family, recorded in the output
adjust_multiplicity <- function(results, method = "holm",
                                family_label = "all reported tests") {
  if (!"p_raw" %in% names(results)) stop("results must contain a p_raw column")
  if (!method %in% stats::p.adjust.methods) {
    stop("method must be one of: ", paste(stats::p.adjust.methods, collapse = ", "))
  }
  testable <- !is.na(results$p_raw)
  results$p_adjusted <- NA_real_
  results$p_adjusted[testable] <- stats::p.adjust(results$p_raw[testable], method = method)
  results$adjustment_method <- method
  results$adjustment_family <- family_label
  results$adjustment_family_size <- sum(testable)
  results
}

#' Run one model per response and collect the group tests.
#'
#' Responses whose design is unusable are reported with the reason rather than
#' dropped, so a missing row in the results is never silent.
#'
#' @param data one row per observation
#' @param responses the measurement columns to test
#' @param ... passed to experiment_design()
#' @param adjust_method multiplicity adjustment method
analyze_responses <- function(data, responses, ..., adjust_method = "holm") {
  rows <- lapply(responses, function(response) {
    design <- tryCatch(
      experiment_design(data, response = response, ...),
      error = function(e) e
    )
    if (inherits(design, "error")) {
      return(data.frame(
        response = response, term = NA_character_, model = NA_character_,
        statistic = NA_real_, df = NA_real_, p_raw = NA_real_,
        n_observations = NA_integer_, n_animals = NA_integer_,
        note = conditionMessage(design), stringsAsFactors = FALSE
      ))
    }
    model <- tryCatch(fit_group_model(design), error = function(e) e)
    if (inherits(model, "error")) {
      validation <- validate_experiment_design(design)
      return(data.frame(
        response = response, term = design$group, model = NA_character_,
        statistic = NA_real_, df = NA_real_, p_raw = NA_real_,
        n_observations = validation$n_observations,
        n_animals = validation$n_animals,
        note = conditionMessage(model), stringsAsFactors = FALSE
      ))
    }
    test_group_effect(model)
  })
  results <- do.call(rbind, rows)
  adjust_multiplicity(
    results, method = adjust_method,
    family_label = sprintf("%d behavioral response(s) in one analysis", length(responses))
  )
}

#' Summarise a design and its chosen model without fitting anything.
#'
#' Useful as a pre-registration artefact: it states what will be fitted and
#' why, before any p-value exists.
describe_analysis_plan <- function(design) {
  validation <- validate_experiment_design(design)
  structure_spec <- model_structure(design)
  list(
    response = design$response,
    group = design$group,
    animal_id = design$animal_id,
    n_observations = validation$n_observations,
    n_animals = validation$n_animals,
    observations_per_animal = c(
      min = validation$min_observations_per_animal,
      max = validation$max_observations_per_animal
    ),
    model_type = structure_spec$type,
    formula = deparse(structure_spec$formula),
    rationale = structure_spec$rationale,
    problems = validation$problems,
    warnings = validation$warnings
  )
}
