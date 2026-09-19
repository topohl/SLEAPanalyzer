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
  # OFT is deliberately absent: it does not use the common assay schema. It is
  # covered by its own test below, against the surface it actually consumes.
  assays <- list(
    NOR = "nor.example.yaml",
    SocP = "socp.example.yaml",
    EPM = "epm.example.yaml"
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

# Read the `config <- list(...)` block out of an assay script without running
# the script, so a test can see the keys it actually accepts.
script_default_config <- function(script_path) {
  src <- readLines(script_path, warn = FALSE)
  start <- grep("^config <- list\\(", src)[1]
  if (is.na(start)) stop("no 'config <- list(' block in ", basename(script_path))
  depth <- 0L
  end <- NA_integer_
  for (i in seq(start, length(src))) {
    opens <- lengths(regmatches(src[i], gregexpr("(", src[i], fixed = TRUE)))
    closes <- lengths(regmatches(src[i], gregexpr(")", src[i], fixed = TRUE)))
    depth <- depth + opens - closes
    if (depth == 0L) { end <- i; break }
  }
  if (is.na(end)) stop("unterminated config block in ", basename(script_path))
  eval(parse(text = paste(src[seq(start, end)], collapse = "\n")))
}

test_that("the OFT example only sets keys the OFT script accepts", {
  skip_if_no_yaml()
  # DLCA_OFT overlays YAML onto its own internal list via apply_config_overlay(),
  # which aborts on any key that list does not define. An example written to
  # the common assay schema therefore parses fine and still cannot be run --
  # which is exactly what shipped before. Check the real contract.
  script <- file.path(repo_root, "02_SLEAPanalzyer", "DLCA_OFT v1.2.0.R")
  skip_if_not(file.exists(script))
  defaults <- script_default_config(script)
  example <- yaml::yaml.load_file(file.path(repo_root, "config", "oft.example.yaml"))

  unknown <- setdiff(names(example), names(defaults))
  expect_equal(
    unknown, character(0),
    info = paste("keys the OFT script would reject:", paste(unknown, collapse = ", "))
  )

  merged <- utils::modifyList(defaults, example)
  expect_true(is.numeric(merged$fps) && merged$fps > 0)
  expect_equal(length(merged$corner_points), 4L)
  expect_true(is.numeric(merged$arena_size_cm) && merged$arena_size_cm > 0)
  # Per-frame displacement, not a speed: a large value silently disables the
  # jump filter, so the shipped example must not carry one.
  expect_lt(merged$max_jump_cm, 50)
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

test_that("apply_config_overlay overrides only the keys it names", {
  skip_if_no_yaml()
  base <- list(fps = 30, input_dir = "orig", threshold = 5, flag = TRUE)
  path <- write_config(c("fps: 60", "threshold: 9"))
  merged <- apply_config_overlay(base, path)
  expect_equal(merged$fps, 60)
  expect_equal(merged$threshold, 9)
  expect_equal(merged$input_dir, "orig")
  expect_true(merged$flag)
})

test_that("apply_config_overlay rejects unknown keys and type changes", {
  skip_if_no_yaml()
  base <- list(fps = 30, label = "a")
  expect_error(
    apply_config_overlay(base, write_config("fpss: 60")),
    "does not use"
  )
  expect_error(
    apply_config_overlay(base, write_config("fps: not-a-number")),
    "changes the type"
  )
})

test_that("apply_config_overlay is a no-op without a configuration", {
  old <- Sys.getenv("SLEAP_ANALYZER_CONFIG", unset = NA_character_)
  Sys.unsetenv("SLEAP_ANALYZER_CONFIG")
  on.exit({
    if (!is.na(old)) Sys.setenv(SLEAP_ANALYZER_CONFIG = old)
  }, add = TRUE)

  base <- list(fps = 30)
  expect_equal(apply_config_overlay(base), base)
  expect_error(apply_config_overlay(base, required = TRUE), "SLEAP_ANALYZER_CONFIG")
})

test_that("the OFT script no longer contains an absolute functions path", {
  script <- readLines(
    file.path(repo_root, "02_SLEAPanalzyer", "DLCA_OFT v1.2.0.R"), warn = FALSE
  )
  functions_line <- grep("functions_file *=", script, value = TRUE)[1]
  expect_false(grepl("[A-Za-z]:/", functions_line))
  expect_true(grepl("DLCAnalyzer_Functions_final.R", functions_line))
})
