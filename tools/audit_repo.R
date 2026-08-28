#!/usr/bin/env Rscript
# Static repository audit, run in CI as a ratchet.
#
# Each pattern below marks a class of problem this repository is actively
# retiring. The audit does not demand zero occurrences today; it records the
# current count as a budget and fails if the count *increases*. That prevents
# a fixed problem from creeping back while migration is still in progress.
#
# When you fix occurrences, lower the corresponding budget in the same commit.

repo_root <- normalizePath(".")

r_files <- function() {
  files <- list.files(repo_root, pattern = "\\.[rR]$", recursive = TRUE, full.names = TRUE)
  files[!grepl("(^|/)(99_deprecated|90_Testing|renv|tools|tests)/", files)]
}

# Blank out double-quoted string literals so that a pattern describing *code*
# is not matched inside a message telling the user what to do.
strip_string_literals <- function(lines) {
  gsub('"[^"]*"', '""', lines)
}

count_pattern <- function(pattern, files, exclude_comments = TRUE,
                          strip_strings = FALSE, perl = FALSE) {
  hits <- list()
  for (path in files) {
    lines <- readLines(path, warn = FALSE)
    if (exclude_comments) lines[grepl("^\\s*#", lines)] <- ""
    if (strip_strings) lines <- strip_string_literals(lines)
    matched <- grep(pattern, lines, perl = perl)
    if (length(matched) > 0) {
      hits[[length(hits) + 1L]] <- data.frame(
        file = sub(paste0("^", repo_root, "/?"), "", path),
        line = matched,
        text = trimws(lines[matched]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(hits) == 0) {
    return(data.frame(file = character(), line = integer(), text = character()))
  }
  do.call(rbind, hits)
}

CHECKS <- list(
  list(
    id = "absolute_paths",
    description = "machine-specific absolute paths in production scripts",
    pattern = "\"[A-Za-z]:/|\"//|\"\\\\\\\\\\\\\\\\",
    budget = 14
  ),
  list(
    # A call is acceptable when it passes an explicit maxgap.
    id = "unbounded_interpolation",
    description = "interpolation without an explicit gap limit",
    pattern = "(na_interpolation|na[.]locf|na[.]approx)[(](?![^)]*maxgap)",
    perl = TRUE,
    strip_strings = TRUE,
    budget = 0
  ),
  list(
    id = "legacy_transition_count",
    description = "CalculateTransitions() used instead of count_entries()",
    pattern = "CalculateTransitions\\(",
    budget = 3
  ),
  list(
    id = "boundary_excluding_zone_test",
    description = "point.in.polygon(...) == 1, which drops boundary points",
    pattern = "point\\.in\\.polygon\\(.*\\)\\s*==\\s*1",
    budget = 0
  ),
  list(
    id = "install_at_runtime",
    description = "package installation during analysis",
    pattern = "install[.]packages[(]",
    strip_strings = TRUE,
    budget = 0
  )
)

failures <- 0L
cat("Static repository audit\n")
cat(strrep("-", 70), "\n")
files <- r_files()

for (check in CHECKS) {
  hits <- count_pattern(
    check$pattern, files,
    strip_strings = isTRUE(check$strip_strings),
    perl = isTRUE(check$perl)
  )
  n <- nrow(hits)
  status <- if (n > check$budget) {
    failures <- failures + 1L
    "OVER BUDGET"
  } else if (n < check$budget) {
    "improved"
  } else {
    "ok"
  }
  cat(sprintf(
    "%-32s %3d / budget %3d  %s\n", check$id, n, check$budget, status
  ))
  if (n > check$budget || n < check$budget) {
    if (n > 0) {
      for (i in seq_len(min(n, 20L))) {
        cat(sprintf("    %s:%d  %s\n", hits$file[i], hits$line[i],
                    substr(hits$text[i], 1, 90)))
      }
      if (n > 20L) cat(sprintf("    ... and %d more\n", n - 20L))
    }
    if (n < check$budget) {
      cat(sprintf(
        "    Budget for '%s' can be lowered to %d.\n", check$id, n
      ))
    }
  }
}

cat(strrep("-", 70), "\n")
if (failures > 0L) {
  cat(sprintf("%d check(s) over budget\n", failures))
  quit(status = 1L)
}
cat("All checks within budget\n")
