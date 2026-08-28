skip_if_no_yaml <- function() {
  testthat::skip_if_not_installed("yaml")
}

write_config <- function(text) {
  path <- tempfile(fileext = ".yaml")
  writeLines(text, path)
  path
}

test_that("merge_config merges nested lists key by key", {
  defaults <- list(fps = 30, qc = list(min_valid = 0.8, max_gap = 5))
  overrides <- list(qc = list(max_gap = 2))
  merged <- merge_config(defaults, overrides)
  expect_equal(merged$fps, 30)
  expect_equal(merged$qc$min_valid, 0.8)
  expect_equal(merged$qc$max_gap, 2)
})

test_that("validate_config reports every problem at once", {
  schema <- list(
    fps = config_field("number", positive = TRUE),
    name = config_field("string"),
    mode = config_field("string", choices = c("a", "b")),
    flag = config_field("logical")
  )
  bad <- list(fps = -1, name = "", mode = "c", flag = "yes")
  err <- tryCatch(validate_config(bad, schema), error = function(e) conditionMessage(e))
  expect_match(err, "fps")
  expect_match(err, "name")
  expect_match(err, "mode")
  expect_match(err, "flag")
})

test_that("validate_config rejects missing required fields and unknown keys", {
  schema <- list(fps = config_field("number", positive = TRUE))
  expect_error(validate_config(list(), schema), "'fps' is required")
  expect_error(validate_config(list(fps = 30, typo = 1), schema), "unknown configuration key")
})

test_that("validate_config enforces vector lengths", {
  schema <- list(range = config_field("number_vector", length = 2))
  expect_error(validate_config(list(range = c(1, 2, 3)), schema), "must have length 2")
  expect_silent(validate_config(list(range = c(1, 2)), schema))
})

test_that("optional fields may be absent", {
  schema <- list(
    fps = config_field("number", positive = TRUE),
    cutoff = config_field("number", required = FALSE)
  )
  expect_silent(validate_config(list(fps = 30), schema))
})

test_that("relative paths resolve against the config file directory", {
  skip_if_no_yaml()
  path <- write_config(c("fps: 30", "input_dir: data/in"))
  config <- load_analysis_config(
    path,
    defaults = list(),
    schema = list(fps = config_field("number", positive = TRUE),
                  input_dir = config_field("path")),
    path_fields = "input_dir"
  )
  expect_equal(config$input_dir, file.path(dirname(normalizePath(path, winslash = "/")), "data/in"))
  expect_true(nzchar(config$.config_path))
})

test_that("absolute paths are left untouched", {
  skip_if_no_yaml()
  absolute <- if (.Platform$OS.type == "windows") "C:/data/in" else "/data/in"
  path <- write_config(c("fps: 30", paste0("input_dir: ", absolute)))
  config <- load_analysis_config(
    path,
    schema = list(fps = config_field("number", positive = TRUE),
                  input_dir = config_field("path")),
    path_fields = "input_dir"
  )
  expect_equal(config$input_dir, absolute)
})

test_that("a missing configuration file is reported clearly", {
  skip_if_no_yaml()
  expect_error(
    load_analysis_config(file.path(tempdir(), "absent.yaml")),
    "Configuration file not found"
  )
})

test_that("the bundled example configurations load and validate", {
  skip_if_no_yaml()
  assays <- list(
    NOR = "nor.example.yaml",
    SocP = "socp.example.yaml",
    EPM = "epm.example.yaml",
    OFT = "oft.example.yaml"
  )
  for (assay in names(assays)) {
    path <- file.path(repo_root, "config", assays[[assay]])
    expect_true(file.exists(path), info = assay)
    config <- load_assay_config(path, assay)
    expect_true(is.numeric(config$fps), info = assay)
    expect_gt(config$fps, 0)
    expect_equal(length(config$arena_corner_names), 4L, info = assay)
  }
})

test_that("the NOR example does not default to the biased legacy detector", {
  skip_if_no_yaml()
  config <- load_assay_config(file.path(repo_root, "config", "nor.example.yaml"), "NOR")
  expect_equal(config$contact_geometry, "radial")
  expect_false(identical(config$contact_geometry, "legacy_asymmetric"))
})

test_that("an invalid assay configuration is rejected", {
  skip_if_no_yaml()
  path <- write_config(c(
    "input_dir: in", "output_dir: out", "fps: -5",
    "arena_width_cm: 49", "arena_height_cm: 49",
    "novel_location_file: novelLoc.txt",
    "contact_geometry: teleport",
    "contact_distance_cm: 4", "body_exclusion_distance_cm: 1",
    "contact_angle_deg: [70, 290]", "proximity_range_cm: [4, 8]",
    "proximity_angle_deg: [90, 270]"
  ))
  err <- tryCatch(load_assay_config(path, "NOR"), error = function(e) conditionMessage(e))
  expect_match(err, "fps")
  expect_match(err, "contact_geometry")
})

test_that("resolve_config_path prefers the environment variable", {
  old <- Sys.getenv("SLEAP_ANALYZER_CONFIG", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("SLEAP_ANALYZER_CONFIG") else
      Sys.setenv(SLEAP_ANALYZER_CONFIG = old)
  }, add = TRUE)

  Sys.setenv(SLEAP_ANALYZER_CONFIG = "/from/env.yaml")
  expect_equal(resolve_config_path(NULL, "NOR"), "/from/env.yaml")
  # An explicit argument still wins.
  expect_equal(resolve_config_path("/explicit.yaml", "NOR"), "/explicit.yaml")
})
