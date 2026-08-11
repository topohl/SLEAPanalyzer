if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("The regression suite requires the testthat package. Install it with install.packages('testthat').")
}

repo_root <- normalizePath(getwd())
if (!file.exists(file.path(repo_root, "02_SLEAPanalzyer", "DLCAnalyzer_Functions_final.R"))) {
  stop("Run tests/testthat.R from the repository root")
}
Sys.setenv(SLEAP_ANALYZER_REPO_ROOT = repo_root)

testthat::test_dir(
  file.path(repo_root, "tests", "testthat"),
  reporter = "summary",
  stop_on_failure = TRUE
)
