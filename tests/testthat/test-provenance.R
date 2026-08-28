test_that("a run manifest records code, environment, config and inputs", {
  input <- tempfile(fileext = ".csv")
  writeLines(c("a,b", "1,2"), input)

  manifest <- run_manifest(
    config = list(fps = 30, contact_distance_cm = 4, .config_dir = "/should/be/dropped"),
    inputs = input,
    packages = c("stats"),
    repo_dir = repo_root
  )

  expect_equal(manifest$schema$output_schema_version, OUTPUT_SCHEMA_VERSION)
  expect_match(manifest$run$timestamp_utc, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$")
  expect_equal(manifest$config$fps, 30)
  # Internal bookkeeping keys must not leak into the manifest.
  expect_null(manifest$config$.config_dir)

  expect_length(manifest$inputs, 1L)
  expect_true(manifest$inputs[[1]]$exists)
  expect_equal(nchar(manifest$inputs[[1]]$md5), 32L)
  expect_gt(manifest$inputs[[1]]$bytes, 0)

  expect_true(nzchar(manifest$environment$r_version))
  expect_true(!is.null(manifest$environment$packages$stats))
})

test_that("the commit of this checkout is recorded", {
  commit <- repository_commit(repo_root)
  expect_true(is.na(commit) || grepl("^[0-9a-f]{40}$", commit))
})

test_that("a missing input file is recorded rather than silently skipped", {
  manifest <- run_manifest(inputs = file.path(tempdir(), "absent.csv"))
  expect_length(manifest$inputs, 1L)
  expect_false(manifest$inputs[[1]]$exists)
})

test_that("the input hash changes when the file changes", {
  input <- tempfile(fileext = ".csv")
  writeLines("original", input)
  first <- run_manifest(inputs = input)$inputs[[1]]$md5
  writeLines("edited", input)
  second <- run_manifest(inputs = input)$inputs[[1]]$md5
  expect_false(identical(first, second))
})

test_that("a manifest round-trips through YAML", {
  testthat::skip_if_not_installed("yaml")
  path <- file.path(tempdir(), "manifest", "run.yaml")
  manifest <- run_manifest(config = list(fps = 30), packages = "stats")
  write_run_manifest(manifest, path)
  expect_true(file.exists(path))

  restored <- yaml::yaml.load_file(path)
  expect_equal(restored$config$fps, 30)
  expect_equal(restored$schema$output_schema_version, OUTPUT_SCHEMA_VERSION)
})

test_that("an unreproducible run is flagged", {
  dirty <- list(code = list(commit = "abc", dirty_working_tree = TRUE))
  expect_warning(assert_reproducible_run(dirty), "uncommitted changes")
  expect_error(assert_reproducible_run(dirty, strict = TRUE), "uncommitted changes")

  untracked <- list(code = list(commit = NA_character_, dirty_working_tree = FALSE))
  expect_warning(assert_reproducible_run(untracked), "not running from a git checkout")

  clean <- list(code = list(commit = "abc", dirty_working_tree = FALSE))
  expect_true(assert_reproducible_run(clean))
})

test_that("package versions are recorded, with NA for absent packages", {
  versions <- package_versions(c("stats", "definitely.not.a.package"))
  expect_true(nzchar(versions$stats))
  expect_true(is.na(versions$definitely.not.a.package))
})
