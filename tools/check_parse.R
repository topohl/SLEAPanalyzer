#!/usr/bin/env Rscript
# Parse every R source file in the repository.
#
# A syntax error in a production script is invisible to the test suite until
# somebody runs that script against real data, which is exactly when it costs
# the most. This check is cheap and catches it immediately.

repo_root <- normalizePath(".")
targets <- list.files(
  repo_root,
  pattern = "\\.[rR]$",
  recursive = TRUE,
  full.names = TRUE
)
# 99_deprecated is retained for provenance only and is not maintained.
targets <- targets[!grepl("(^|/)(99_deprecated|renv)/", targets)]

failures <- character()
for (path in targets) {
  result <- tryCatch({
    parse(path)
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(result)) {
    failures <- c(failures, paste0(
      sub(paste0("^", repo_root, "/?"), "", path), ": ", result
    ))
  }
}

cat(sprintf("Parsed %d R files\n", length(targets)))
if (length(failures) > 0) {
  cat("\nParse failures:\n")
  cat(paste0("  - ", failures, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("All R files parse cleanly\n")
