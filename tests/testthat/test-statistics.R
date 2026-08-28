# The statistical layer: design validation and model selection.

source_statistics <- function() {
  source(file.path(repo_root, "03_statistics", "design.R"), local = parent.frame())
  source(file.path(repo_root, "03_statistics", "models.R"), local = parent.frame())
}
source_statistics()

one_session_per_animal <- function(n = 20) {
  data.frame(
    ID = paste0("m", seq_len(n)),
    group = rep(c("control", "stress"), length.out = n),
    sex = rep(c("F", "M"), each = n / 2),
    batch = rep(c("B1", "B2"), length.out = n),
    center_time = c(rnorm(n / 2, 30, 5), rnorm(n / 2, 20, 5)),
    stringsAsFactors = FALSE
  )
}

repeated_sessions <- function(n_animals = 10, sessions = 3) {
  data.frame(
    ID = rep(paste0("m", seq_len(n_animals)), each = sessions),
    group = rep(rep(c("control", "stress"), length.out = n_animals), each = sessions),
    phase = rep(paste0("S", seq_len(sessions)), times = n_animals),
    center_time = rnorm(n_animals * sessions, 25, 5),
    stringsAsFactors = FALSE
  )
}

test_that("one observation per animal yields a fixed-effects model", {
  set.seed(1)
  design <- experiment_design(
    one_session_per_animal(), response = "center_time",
    group = "group", animal_id = "ID"
  )
  structure_spec <- model_structure(design)

  # This is the case the old extension script got wrong: with one session per
  # animal, lme4 refuses (1|ID) outright because the number of groups equals
  # the number of observations.
  expect_false(structure_spec$repeated_measures)
  expect_equal(structure_spec$type, "linear model")
  expect_length(structure_spec$random_effects, 0)
  expect_match(structure_spec$rationale, "unidentifiable")
})

test_that("repeated observations per animal yield a mixed model", {
  set.seed(2)
  design <- experiment_design(
    repeated_sessions(), response = "center_time",
    group = "group", animal_id = "ID", within = "phase"
  )
  structure_spec <- model_structure(design)
  expect_true(structure_spec$repeated_measures)
  expect_equal(structure_spec$type, "linear mixed model")
  expect_equal(structure_spec$random_effects, "(1|ID)")
  expect_true("phase" %in% structure_spec$fixed_effects)
})

test_that("the model structure does not depend on installed packages", {
  set.seed(3)
  design <- experiment_design(
    one_session_per_animal(), response = "center_time",
    group = "group", animal_id = "ID"
  )
  # model_structure() consults only the data. Calling it twice with the same
  # data must give the same answer regardless of the library state.
  # Formulas carry their environment, so compare the parts that matter.
  comparable <- function(spec) {
    spec$formula <- deparse(spec$formula)
    spec
  }
  expect_identical(comparable(model_structure(design)), comparable(model_structure(design)))
  expect_false(model_structure(design)$repeated_measures)

  # The decisive property: the choice is a function of the data alone, so it
  # cannot change because lmerTest was installed or removed.
  expect_equal(
    model_structure(design)$repeated_measures,
    max(observations_per_animal(design)) > 1
  )
})

test_that("a missing mixed-model package is an error, not a silent fallback", {
  set.seed(4)
  design <- experiment_design(
    repeated_sessions(), response = "center_time",
    group = "group", animal_id = "ID", within = "phase"
  )
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    model <- fit_group_model(design)
    expect_equal(model$structure$type, "linear mixed model")
  } else {
    expect_error(fit_group_model(design), "lmerTest")
  }
})

test_that("a design with one group level is rejected", {
  data <- one_session_per_animal()
  data$group <- "control"
  design <- experiment_design(data, response = "center_time",
                              group = "group", animal_id = "ID")
  validation <- validate_experiment_design(design)
  expect_false(validation$ok)
  expect_match(validation$problems, "fewer than two levels", all = FALSE)
  expect_error(fit_group_model(design), "does not support an analysis")
})

test_that("complete confounding between batch and group is reported", {
  data <- one_session_per_animal(20)
  # Every animal in B1 is a control and every animal in B2 is stressed.
  data$batch <- ifelse(data$group == "control", "B1", "B2")
  design <- experiment_design(
    data, response = "center_time", group = "group",
    animal_id = "ID", batch = "batch"
  )
  validation <- validate_experiment_design(design)
  expect_false(validation$ok)
  expect_match(validation$problems, "completely confounded", all = FALSE)
})

test_that("an animal appearing in two groups without a within factor is rejected", {
  data <- repeated_sessions(6, 2)
  data$group <- rep(c("control", "stress"), times = 6)
  design <- experiment_design(data, response = "center_time",
                              group = "group", animal_id = "ID")
  validation <- validate_experiment_design(design)
  expect_false(validation$ok)
  expect_match(validation$problems, "more than one level", all = FALSE)
})

test_that("undeclared repeated measures produce a warning, not silence", {
  set.seed(5)
  data <- repeated_sessions(8, 3)
  design <- experiment_design(data, response = "center_time",
                              group = "group", animal_id = "ID")
  validation <- validate_experiment_design(design)
  expect_true(validation$ok)
  expect_match(validation$warnings, "within-animal factor", all = FALSE)
})

test_that("time bins are not treated as independent animals", {
  set.seed(6)
  # Ten animals, six one-minute bins each. Treating the 60 rows as independent
  # would quadruple the apparent sample size.
  data <- repeated_sessions(10, 6)
  names(data)[names(data) == "phase"] <- "time_bin"
  design <- experiment_design(
    data, response = "center_time", group = "group",
    animal_id = "ID", within = "time_bin"
  )
  validation <- validate_experiment_design(design)
  expect_equal(validation$n_observations, 60L)
  expect_equal(validation$n_animals, 10L)
  # The reported sample size distinguishes animals from observations.
  expect_true(model_structure(design)$repeated_measures)
})

test_that("multiplicity adjustment records the family and its size", {
  results <- data.frame(
    response = c("a", "b", "c"),
    p_raw = c(0.01, 0.04, 0.20),
    stringsAsFactors = FALSE
  )
  adjusted <- adjust_multiplicity(results, method = "holm", family_label = "three metrics")
  expect_equal(adjusted$p_adjusted, p.adjust(results$p_raw, method = "holm"))
  expect_equal(unique(adjusted$adjustment_method), "holm")
  expect_equal(unique(adjusted$adjustment_family), "three metrics")
  expect_equal(unique(adjusted$adjustment_family_size), 3L)
})

test_that("failed models do not inflate the multiplicity family", {
  results <- data.frame(
    response = c("a", "b", "c", "d"),
    p_raw = c(0.01, NA, 0.04, 0.03),
    stringsAsFactors = FALSE
  )
  adjusted <- adjust_multiplicity(results)
  # Three tests were actually performed, so the family size is three.
  expect_equal(unique(adjusted$adjustment_family_size), 3L)
  expect_equal(
    adjusted$p_adjusted[!is.na(adjusted$p_raw)],
    p.adjust(c(0.01, 0.04, 0.03), method = "holm")
  )
  expect_true(is.na(adjusted$p_adjusted[2]))
})

test_that("adjust_multiplicity validates its arguments", {
  results <- data.frame(p_raw = c(0.1, 0.2))
  expect_error(adjust_multiplicity(data.frame(x = 1)), "p_raw")
  expect_error(adjust_multiplicity(results, method = "magic"), "method must be one of")
})

test_that("analyze_responses reports unusable responses instead of dropping them", {
  set.seed(7)
  data <- one_session_per_animal()
  data$all_missing <- NA_real_
  results <- analyze_responses(
    data, responses = c("center_time", "all_missing"),
    group = "group", animal_id = "ID"
  )
  expect_equal(nrow(results), 2L)
  expect_true(is.na(results$p_raw[results$response == "all_missing"]))
  expect_true(nzchar(results$note[results$response == "all_missing"]))
  expect_false(is.na(results$p_raw[results$response == "center_time"]))
})

test_that("a group difference is detected and reported with sample sizes", {
  set.seed(8)
  data <- one_session_per_animal(40)
  data$center_time <- ifelse(data$group == "control", 30, 15) + rnorm(40, 0, 2)
  results <- analyze_responses(
    data, responses = "center_time", group = "group", animal_id = "ID"
  )
  expect_equal(results$model, "linear model")
  expect_lt(results$p_raw, 0.001)
  expect_equal(results$n_observations, 40L)
  expect_equal(results$n_animals, 40L)
})

test_that("describe_analysis_plan states the model before any p-value exists", {
  set.seed(9)
  design <- experiment_design(
    one_session_per_animal(), response = "center_time",
    group = "group", animal_id = "ID"
  )
  plan <- describe_analysis_plan(design)
  expect_equal(plan$model_type, "linear model")
  expect_match(plan$formula, "center_time ~ group")
  expect_match(plan$rationale, "exactly one observation")
  expect_equal(plan$observations_per_animal[["max"]], 1L)
})

test_that("a non-numeric response is rejected", {
  data <- one_session_per_animal()
  data$label <- "x"
  expect_error(
    experiment_design(data, response = "label", group = "group", animal_id = "ID"),
    "must be numeric"
  )
})
