#!/usr/bin/env Rscript
# Compare two sets of scoring outputs, table by table and cell by cell.
#
# Two modes:
#   tree     every CSV file under OLD and under NEW, matched by its path relative
#            to each folder (an output folder before and after a rescoring, or the
#            outputs of two scoring runs).
#   release  the CSV files of a release's assay_summaries/ against the live files
#            the release build copies them from (tools/build_exp9_behavior_bundle.ps1):
#            EPM, NOR, NOR_HAB and SocP from the assay output root, OFT from
#            <stage root>/OFT/<B>/output_v1.2.0. Old is the release, new the live file.
#            A live file on the build's list (release_build_files()) that the
#            release lacks is listed as only in new.
#
# Usage:
#   Rscript tools/compare_outputs.R OLD NEW [--out diffs.csv] [--allow allow.csv]
#   Rscript tools/compare_outputs.R --release RELEASE [--output-root DIR] [--stage-root DIR]
#                                   [--out diffs.csv] [--allow allow.csv]
#   More options: --tolerance 1e-9, --key "combined_bins.csv=Code,bin_s" (repeatable),
#   --prefix NOR/B1 (tree mode: put this in front of every file path).
#   RELEASE is the release folder or its assay_summaries/. --output-root and
#   --stage-root default to EXP9_SLEAP_ASSAY_OUTPUT_ROOT and EXP9_SLEAP_STAGE_ROOT,
#   else to the release build's defaults.
#   Or source() this file and call compare_output_trees() or compare_release().
#
# Rows are matched by key: the columns named for that file (--key, or the keys of
# the multi-row tables in compare_default_keys()), else the first of Code and file
# that is unique in both tables, else Code or file, else the first column. Where a
# key repeats in a table, every row also gets its occurrence number ("A1L9 #2").
# An empty table (write.csv(data.frame()) writes one line, "") has no columns; its
# rows count as removed or added under the other table's key. Two cells are equal
# when their text is equal, or when both are numbers within the tolerance relative
# to the larger of 1 and their size.
#
# Each file is reported as identical, bytes differ but every cell is equal,
# changed (columns added or removed, rows added or removed, changed cells), or
# only in old or only in new. --out writes every difference as one row of a CSV:
# file, change, key, column, old, new, expected.
#
# Allowlist (--allow): a CSV with the columns file, key and column (others, such
# as a reason, are ignored). file and column are glob patterns, empty for any;
# a file pattern without "/" matches the file name alone. key is a list separated
# by ";", empty for any; each entry is a glob, or a code that also matches the
# keys that start with it (A1L9 matches A1L9_NOV and "A1L9 | 60"). A row allows a
# difference when all three match. A row that names keys does not allow an added
# or removed column, and a row that names a column does not allow an added or
# removed row or file. Every difference no row allows is unexpected.
#
# Exit status: 0 done, 1 differences outside the allowlist, 2 usage or input error.

# Keys of the output tables with several rows per animal, by file name (glob).
compare_default_keys <- function() {
  list(
    "combined_bins.csv" = c("Code", "bin_s"),
    "combined_first_last_epoch_metrics.csv" = c("Code", "epoch"),
    "combined_qc.csv" = c("file", "stage", "point"),
    "landmark_repairs.csv" = c("file", "landmark"),
    "landmark_repair_models.csv" = c("file", "target", "rank"),
    "on_cup_episodes.csv" = c("file", "object", "start_frame"),
    "landmark_repair_episodes.csv" = c("file", "check", "target", "start_frame"),
    "hidden_episodes.csv" = c("file", "target", "start_frame")
  )
}

# The CSV files the release build copies, by assay: its batches, its phase folders
# and the file names. Keep in step with tools/build_exp9_behavior_bundle.ps1
# (lines 155-186), which also copies each run_manifest.yaml.
release_build_files <- function() {
  batches <- paste0("B", 1:6)
  list(
    EPM = list(batches = batches, files = c("Report.csv", "tracking_qc.csv")),
    NOR = list(batches = batches, files = c("combined_output.csv", "tracking_qc.csv")),
    # Batch 1's habituation trials were never tracked.
    NOR_HAB = list(batches = batches[-1], files = c("combined_output.csv", "tracking_qc.csv")),
    SocP = list(batches = batches, phases = c("HAB", "S1", "S2"),
                files = c("combined_output.csv", "tracking_qc.csv", "on_cup_episodes.csv", "landmark_repairs.csv",
                          "landmark_repair_episodes.csv", "landmark_repair_models.csv")),
    OFT = list(batches = batches, files = c("combined_summary.csv", "combined_qc.csv", "combined_bins.csv",
                                            "combined_enhanced_metrics.csv", "combined_first_last_epoch_metrics.csv"))
  )
}

# The release build's default roots, used when neither an option nor the
# environment variable names them.
compare_release_defaults <- function() {
  list(
    output_root = "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris/sleap_output_all",
    stage_root = "C:/Users/topohl/Documents/exp9_sleap_staging"
  )
}

# ---------------------------------------------------------------------------
# Pairs of files to compare
# ---------------------------------------------------------------------------

list_csv <- function(dir) {
  if (!dir.exists(dir)) return(character())
  sort(list.files(dir, pattern = "\\.[cC][sS][vV]$", recursive = TRUE))
}

# Every CSV file under either folder, matched by relative path.
tree_pairs <- function(old_dir, new_dir, prefix = "") {
  if (!dir.exists(old_dir) && !dir.exists(new_dir)) {
    stop("neither folder exists: ", old_dir, ", ", new_dir)
  }
  a <- list_csv(old_dir)
  b <- list_csv(new_dir)
  rel <- sort(union(a, b))
  data.frame(
    file = if (nzchar(prefix)) paste(sub("/+$", "", prefix), rel, sep = "/") else rel,
    old = ifelse(rel %in% a, file.path(old_dir, rel), NA_character_),
    new = ifelse(rel %in% b, file.path(new_dir, rel), NA_character_),
    stringsAsFactors = FALSE
  )
}

# Every path the release build writes under assay_summaries/ (CSV files only).
release_build_paths <- function() {
  spec <- release_build_files()
  unlist(lapply(names(spec), function(assay) {
    s <- spec[[assay]]
    folders <- if (is.null(s$phases)) s$batches else as.vector(t(outer(s$batches, s$phases, paste, sep = "/")))
    as.vector(t(outer(paste(assay, folders, sep = "/"), s$files, paste, sep = "/")))
  }), use.names = FALSE)
}

# The live file a release file is copied from.
release_live_path <- function(rel, output_root, stage_root) {
  p <- strsplit(rel, "/", fixed = TRUE)[[1]]
  if (p[1] == "OFT") {
    do.call(file.path, as.list(c(stage_root, "OFT", p[2], "output_v1.2.0", p[-(1:2)])))
  } else {
    file.path(output_root, rel)
  }
}

# The release's CSV files and the live files the release build copies them from,
# plus each live file on the build's list that the release lacks.
release_pairs <- function(release, output_root, stage_root) {
  summaries <- if (basename(sub("[/\\\\]+$", "", release)) == "assay_summaries") release else
    file.path(release, "assay_summaries")
  if (!dir.exists(summaries)) stop("no assay_summaries folder: ", summaries)
  rel <- list_csv(summaries)
  if (length(rel) == 0) stop("no CSV files in ", summaries)
  unknown <- setdiff(unique(sub("/.*$", "", rel)), names(release_build_files()))
  if (length(unknown) > 0) {
    stop("the release holds folder(s) the release build does not write: ", paste(unknown, collapse = ", "))
  }
  live_of <- function(files) vapply(files, release_live_path, "", output_root, stage_root, USE.NAMES = FALSE)
  lacking <- setdiff(release_build_paths(), rel)
  files <- sort(c(rel, lacking[file.exists(live_of(lacking))]))
  live <- live_of(files)
  data.frame(file = files, old = ifelse(files %in% rel, file.path(summaries, files), NA_character_),
             new = ifelse(file.exists(live), live, NA_character_), stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# One pair of tables
# ---------------------------------------------------------------------------

# A CSV file with every cell as text, exactly as written ("NA" stays "NA"). A file
# whose lines are all blank or "" (write.csv(data.frame())) is an empty table.
read_csv_cells <- function(path) {
  txt <- gsub("^\ufeff|[ \t\r]", "", readLines(path, warn = FALSE, encoding = "UTF-8"), useBytes = TRUE)
  if (all(txt %in% c("", "\"\""))) return(data.frame())
  x <- utils::read.csv(path, colClasses = "character", check.names = FALSE, na.strings = character(),
                       strip.white = FALSE, comment.char = "", row.names = NULL, fileEncoding = "UTF-8-BOM")
  names(x) <- make.unique(names(x))
  x
}

# The key columns for one file: a given key, else the first of Code and file that is
# unique in both tables, else Code or file, else the first column. When one table
# has no columns (an empty table), the other table's columns are the candidates.
choose_key <- function(file, old, new, keys) {
  cols <- if (!ncol(old)) names(new) else if (!ncol(new)) names(old) else intersect(names(old), names(new))
  for (pattern in names(keys)) {
    target <- if (grepl("/", pattern, fixed = TRUE)) file else basename(file)
    k <- keys[[pattern]]
    if (grepl(utils::glob2rx(pattern), target) && all(k %in% cols)) return(k)
  }
  candidates <- intersect(c("Code", "file"), cols)
  unique_in <- function(df, k) !k %in% names(df) || !anyDuplicated(df[[k]])
  for (k in candidates) if (unique_in(old, k) && unique_in(new, k)) return(k)
  if (length(candidates)) return(candidates[1])
  if (length(cols) == 0) character() else cols[1]
}

row_keys <- function(df, key) {
  if (nrow(df) == 0) return(character())
  if (length(key) == 0) return(as.character(seq_len(nrow(df))))
  if (length(key) == 1) return(df[[key]])
  do.call(paste, c(unname(as.list(df[key])), sep = " | "))
}

with_occurrence <- function(k) {
  if (length(k) == 0) return(k)
  paste0(k, " #", stats::ave(seq_along(k), k, FUN = seq_along))
}

# TRUE where two text cells hold the same value.
cells_equal <- function(x, y, tolerance) {
  same <- x == y
  todo <- which(!same)
  if (length(todo) > 0) {
    a <- suppressWarnings(as.numeric(x[todo]))
    b <- suppressWarnings(as.numeric(y[todo]))
    num <- !is.na(a) & !is.na(b)
    finite <- num & is.finite(a) & is.finite(b)
    close <- rep(FALSE, length(todo))
    close[finite] <- abs(a[finite] - b[finite]) <= tolerance * pmax(1, abs(a[finite]), abs(b[finite]))
    close[num & !finite] <- a[num & !finite] == b[num & !finite]
    same[todo] <- close
  }
  same
}

# Difference rows; their number is the length of the vectors given (one row when none is given).
diff_rows <- function(file, change, key = NULL, column = NULL, old = NULL, new = NULL) {
  given <- Filter(Negate(is.null), list(key, column, old, new))
  n <- if (length(given)) max(lengths(given)) else 1L
  if (n == 0) return(empty_diffs())
  fill <- function(x) if (is.null(x)) rep("", n) else rep(x, length.out = n)
  data.frame(file = rep(file, n), change = rep(change, n), key = fill(key), column = fill(column),
             old = fill(old), new = fill(new), stringsAsFactors = FALSE)
}

empty_diffs <- function() {
  data.frame(file = character(), change = character(), key = character(), column = character(),
             old = character(), new = character(), stringsAsFactors = FALSE)
}

file_row <- function(file, status, rows_old = NA_integer_, rows_new = NA_integer_, key = "",
                     columns_added = 0L, columns_removed = 0L, rows_added = 0L, rows_removed = 0L,
                     rows_changed = 0L, cells_changed = 0L) {
  data.frame(file = file, status = status, rows_old = rows_old, rows_new = rows_new, key = key,
             columns_added = columns_added, columns_removed = columns_removed, rows_added = rows_added,
             rows_removed = rows_removed, rows_changed = rows_changed, cells_changed = cells_changed,
             stringsAsFactors = FALSE)
}

# Compare one file. Returns list(file = one summary row, diffs = its differences).
compare_csv_pair <- function(file, old_path, new_path, tolerance = 1e-9, keys = compare_default_keys()) {
  have_old <- !is.na(old_path) && file.exists(old_path)
  have_new <- !is.na(new_path) && file.exists(new_path)
  if (!have_old && !have_new) stop("neither copy of ", file, " exists")
  read <- function(path) tryCatch(read_csv_cells(path), error = function(e) e)
  if (!have_old || !have_new) {
    present <- read(if (have_new) new_path else old_path)
    n <- if (is.data.frame(present)) nrow(present) else NA_integer_
    status <- if (have_new) "only_new" else "only_old"
    return(list(
      file = file_row(file, status, rows_old = if (have_old) n else NA_integer_,
                      rows_new = if (have_new) n else NA_integer_),
      diffs = diff_rows(file, if (have_new) "file_added" else "file_removed")
    ))
  }
  if (identical(unname(tools::md5sum(old_path)), unname(tools::md5sum(new_path)))) {
    return(list(file = file_row(file, "identical"), diffs = empty_diffs()))
  }
  old <- read(old_path)
  new <- read(new_path)
  if (!is.data.frame(old) || !is.data.frame(new)) {
    problem <- conditionMessage(if (is.data.frame(old)) new else old)
    return(list(file = file_row(file, "unreadable"),
                diffs = diff_rows(file, "unreadable", old = if (is.data.frame(old)) "" else problem,
                                  new = if (is.data.frame(new)) "" else problem)))
  }

  key <- choose_key(file, old, new, keys)
  ko <- row_keys(old, key)
  kn <- row_keys(new, key)
  if (anyDuplicated(ko) || anyDuplicated(kn)) {
    ko <- with_occurrence(ko)
    kn <- with_occurrence(kn)
  }
  key_label <- if (length(key)) paste(key, collapse = "+") else "row number"

  # An empty table has no columns to add or remove; its rows are what changes.
  shaped <- ncol(old) > 0 && ncol(new) > 0
  added_cols <- if (shaped) setdiff(names(new), names(old)) else character()
  removed_cols <- if (shaped) setdiff(names(old), names(new)) else character()
  added_rows <- kn[!kn %in% ko]
  removed_rows <- ko[!ko %in% kn]
  common <- kn[kn %in% ko]
  io <- match(common, ko)
  inew <- match(common, kn)

  cell_parts <- list()
  for (j in setdiff(intersect(names(new), names(old)), key)) {
    x <- old[[j]][io]
    y <- new[[j]][inew]
    bad <- which(!cells_equal(x, y, tolerance))
    if (length(bad) > 0) {
      cell_parts[[length(cell_parts) + 1]] <- data.frame(
        pos = inew[bad], col = match(j, names(new)), key = common[bad], column = j,
        old = x[bad], new = y[bad], stringsAsFactors = FALSE
      )
    }
  }
  cells <- if (length(cell_parts)) do.call(rbind, cell_parts) else NULL
  if (!is.null(cells)) cells <- cells[order(cells$pos, cells$col), , drop = FALSE]

  diffs <- rbind(
    diff_rows(file, "column_added", column = added_cols),
    diff_rows(file, "column_removed", column = removed_cols),
    diff_rows(file, "row_added", key = added_rows),
    diff_rows(file, "row_removed", key = removed_rows),
    if (is.null(cells)) empty_diffs() else
      diff_rows(file, "cell", key = cells$key, column = cells$column, old = cells$old, new = cells$new)
  )
  status <- if (nrow(diffs) == 0) "equal_cells" else "changed"
  list(
    file = file_row(file, status, rows_old = nrow(old), rows_new = nrow(new), key = key_label,
                    columns_added = length(added_cols), columns_removed = length(removed_cols),
                    rows_added = length(added_rows), rows_removed = length(removed_rows),
                    rows_changed = if (is.null(cells)) 0L else length(unique(cells$key)),
                    cells_changed = if (is.null(cells)) 0L else nrow(cells)),
    diffs = diffs
  )
}

# ---------------------------------------------------------------------------
# Allowlist
# ---------------------------------------------------------------------------

read_allowlist <- function(path) {
  if (!file.exists(path)) stop("allowlist not found: ", path)
  allow <- utils::read.csv(path, colClasses = "character", check.names = FALSE, na.strings = character(),
                           fileEncoding = "UTF-8-BOM")
  missing <- setdiff(c("file", "key", "column"), names(allow))
  if (length(missing) > 0) stop("allowlist ", path, " lacks column(s): ", paste(missing, collapse = ", "))
  for (j in c("file", "key", "column")) allow[[j]] <- trimws(allow[[j]])
  allow[, c("file", "key", "column"), drop = FALSE]
}

key_matches <- function(keys, entry) {
  if (grepl("[*?]", entry)) return(grepl(utils::glob2rx(entry), keys))
  rest <- substring(keys, nchar(entry) + 1)
  keys == entry | (startsWith(keys, entry) & grepl("^[^A-Za-z0-9]", rest))
}

# TRUE for each difference some allowlist row allows.
allowed_differences <- function(diffs, allow) {
  ok <- rep(FALSE, nrow(diffs))
  if (nrow(diffs) == 0 || nrow(allow) == 0) return(ok)
  for (i in seq_len(nrow(allow))) {
    rule <- allow[i, ]
    hit <- rep(TRUE, nrow(diffs))
    if (nzchar(rule$file)) {
      target <- if (grepl("/", rule$file, fixed = TRUE)) diffs$file else basename(diffs$file)
      hit <- hit & grepl(utils::glob2rx(rule$file), target)
    }
    entries <- trimws(strsplit(rule$key, ";", fixed = TRUE)[[1]])
    entries <- entries[nzchar(entries)]
    if (length(entries) > 0) {
      hit <- hit & nzchar(diffs$key) & Reduce(`|`, lapply(entries, function(e) key_matches(diffs$key, e)))
    }
    if (nzchar(rule$column)) {
      hit <- hit & nzchar(diffs$column) & grepl(utils::glob2rx(rule$column), diffs$column)
    }
    ok <- ok | hit
  }
  ok
}

# ---------------------------------------------------------------------------
# Whole comparisons
# ---------------------------------------------------------------------------

# Compare every pair (data frame with file, old, new). allow: NULL, a path or a data frame.
compare_pairs <- function(pairs, tolerance = 1e-9, keys = list(), allow = NULL) {
  keys <- c(keys, compare_default_keys()[setdiff(names(compare_default_keys()), names(keys))])
  if (is.character(allow)) allow <- read_allowlist(allow)
  parts <- lapply(seq_len(nrow(pairs)), function(i) {
    compare_csv_pair(pairs$file[i], pairs$old[i], pairs$new[i], tolerance, keys)
  })
  files <- if (length(parts)) do.call(rbind, lapply(parts, `[[`, "file")) else file_row("", "")[0, ]
  diffs <- if (length(parts)) do.call(rbind, lapply(parts, `[[`, "diffs")) else empty_diffs()
  diffs$expected <- if (is.null(allow)) rep(NA, nrow(diffs)) else allowed_differences(diffs, allow)
  files$unexpected <- if (is.null(allow)) rep(NA_integer_, nrow(files)) else
    vapply(files$file, function(f) sum(!diffs$expected[diffs$file == f]), 0L, USE.NAMES = FALSE)
  rownames(files) <- NULL
  rownames(diffs) <- NULL
  list(files = files, diffs = diffs, tolerance = tolerance, allowlist = !is.null(allow),
       unexpected = if (is.null(allow)) 0L else sum(!diffs$expected))
}

#' Compare every CSV file under two folders.
compare_output_trees <- function(old_dir, new_dir, prefix = "", tolerance = 1e-9, keys = list(), allow = NULL) {
  compare_pairs(tree_pairs(old_dir, new_dir, prefix), tolerance, keys, allow)
}

#' Compare a release's assay_summaries with the live files the release build copies.
compare_release <- function(release, output_root, stage_root, tolerance = 1e-9, keys = list(), allow = NULL) {
  compare_pairs(release_pairs(release, output_root, stage_root), tolerance, keys, allow)
}

#' Write every difference as one row of a CSV.
write_comparison <- function(result, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result$diffs, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}

#' Print the per-file summary, the changed rows of each changed file and the unexpected differences.
print_comparison <- function(result, max_rows = 12, max_columns = 10, max_unexpected = 30) {
  f <- result$files
  count <- function(s) sum(f$status == s)
  cat(sprintf("Compared %d CSV files\n", nrow(f)))
  cat(sprintf("  identical                 %d\n", count("identical")))
  cat(sprintf("  bytes differ, cells equal %d (within %g)\n", count("equal_cells"), result$tolerance))
  cat(sprintf("  changed                   %d\n", count("changed")))
  cat(sprintf("  only in old               %d\n", count("only_old")))
  cat(sprintf("  only in new               %d\n", count("only_new")))
  if (count("unreadable")) cat(sprintf("  unreadable                %d\n", count("unreadable")))
  d <- result$diffs
  for (i in which(f$status != "identical" & f$status != "equal_cells")) {
    r <- f[i, ]
    if (r$status != "changed") {
      label <- c(only_old = "ONLY IN OLD", only_new = "ONLY IN NEW", unreadable = "UNREADABLE")[[r$status]]
      cat(sprintf("%s %s\n", label, r$file))
      next
    }
    mine <- d[d$file == r$file, , drop = FALSE]
    bits <- sprintf("rows %d -> %d (key %s)", r$rows_old, r$rows_new, r$key)
    for (ch in c("column_added", "column_removed")) {
      cols <- mine$column[mine$change == ch]
      if (length(cols)) bits <- c(bits, sprintf("%s: %s", sub("_", "s ", ch), paste(cols, collapse = ", ")))
    }
    if (r$rows_added) bits <- c(bits, sprintf("rows added: %d", r$rows_added))
    if (r$rows_removed) bits <- c(bits, sprintf("rows removed: %d", r$rows_removed))
    bits <- c(bits, sprintf("changed cells: %d in %d rows", r$cells_changed, r$rows_changed))
    cat(sprintf("CHANGED %s: %s\n", r$file, paste(bits, collapse = "; ")))
    cells <- mine[mine$change == "cell", , drop = FALSE]
    shown <- utils::head(unique(cells$key), max_rows)
    for (k in shown) {
      cols <- cells$column[cells$key == k]
      more <- if (length(cols) > max_columns) sprintf(" ... (%d columns)", length(cols)) else ""
      cat(sprintf("    %s: %s%s\n", k, paste(utils::head(cols, max_columns), collapse = ", "), more))
    }
    if (length(unique(cells$key)) > max_rows) cat(sprintf("    ... %d more rows\n", length(unique(cells$key)) - max_rows))
    rows <- mine[mine$change %in% c("row_added", "row_removed"), , drop = FALSE]
    for (j in utils::head(seq_len(nrow(rows)), max_rows)) {
      cat(sprintf("    %s: %s\n", rows$key[j], sub("_", " ", rows$change[j])))
    }
  }
  if (result$allowlist) {
    bad <- d[!d$expected, , drop = FALSE]
    cat(sprintf("Allowlist: %d differences, %d expected, %d unexpected\n", nrow(d), sum(d$expected), nrow(bad)))
    for (i in utils::head(seq_len(nrow(bad)), max_unexpected)) {
      b <- bad[i, ]
      cat(sprintf("  UNEXPECTED %s %s %s %s %s\n", b$file, sub("_", " ", b$change),
                  if (nzchar(b$key)) b$key else "-", if (nzchar(b$column)) b$column else "-",
                  if (b$change == "cell") paste(b$old, "->", b$new) else ""))
    }
    if (nrow(bad) > max_unexpected) cat(sprintf("  ... %d more\n", nrow(bad) - max_unexpected))
  }
  invisible(result)
}

# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

# Split command-line arguments into positional ones and --name value options.
# Options in `repeatable` collect every value given; `flags` take no value.
parse_tool_args <- function(argv, values = character(), flags = character(), repeatable = character()) {
  out <- list(positional = character())
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (startsWith(a, "--")) {
      name <- sub("^--", "", a)
      value <- NULL
      if (grepl("=", name, fixed = TRUE)) {
        value <- sub("^[^=]*=", "", name)
        name <- sub("=.*$", "", name)
      }
      if (name %in% flags) {
        if (!is.null(value)) stop("--", name, " takes no value")
        out[[name]] <- TRUE
      } else if (name %in% c(values, repeatable)) {
        if (is.null(value)) {
          i <- i + 1
          if (i > length(argv)) stop("--", name, " needs a value")
          value <- argv[i]
        }
        out[[name]] <- if (name %in% repeatable) c(out[[name]], value) else value
      } else {
        stop("unknown option --", name)
      }
    } else {
      out$positional <- c(out$positional, a)
    }
    i <- i + 1
  }
  out
}

# --key "pattern=col1,col2" options as a named list.
parse_key_options <- function(specs) {
  keys <- list()
  for (s in specs) {
    if (!grepl("=", s, fixed = TRUE)) stop("--key needs the form pattern=column[,column]: ", s)
    cols <- trimws(strsplit(sub("^[^=]*=", "", s), ",", fixed = TRUE)[[1]])
    keys[[trimws(sub("=.*$", "", s))]] <- cols[nzchar(cols)]
  }
  keys
}

compare_outputs_usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript tools/compare_outputs.R OLD NEW [--out diffs.csv] [--allow allow.csv] [--prefix P]",
    "  Rscript tools/compare_outputs.R --release RELEASE [--output-root DIR] [--stage-root DIR]",
    "                                  [--out diffs.csv] [--allow allow.csv]",
    "  Options: --tolerance 1e-9, --key \"file.csv=col1,col2\" (repeatable)",
    sep = "\n"), "\n")
}

compare_outputs_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  tryCatch({
    o <- parse_tool_args(argv, values = c("out", "allow", "tolerance", "release", "output-root", "stage-root", "prefix"),
                         flags = "help", repeatable = "key")
    if (isTRUE(o$help)) {
      compare_outputs_usage()
      return(0L)
    }
    tolerance <- if (is.null(o$tolerance)) 1e-9 else as.numeric(o$tolerance)
    if (!is.finite(tolerance) || tolerance < 0) stop("--tolerance must be a number >= 0")
    keys <- parse_key_options(o$key)
    allow <- if (is.null(o$allow)) NULL else read_allowlist(o$allow)
    if (!is.null(o$release)) {
      if (length(o$positional) > 0) stop("--release takes no OLD and NEW folders")
      defaults <- compare_release_defaults()
      from_env <- function(name, default) if (nzchar(Sys.getenv(name))) Sys.getenv(name) else default
      output_root <- o[["output-root"]]
      if (is.null(output_root)) output_root <- from_env("EXP9_SLEAP_ASSAY_OUTPUT_ROOT", defaults$output_root)
      stage_root <- o[["stage-root"]]
      if (is.null(stage_root)) stage_root <- from_env("EXP9_SLEAP_STAGE_ROOT", defaults$stage_root)
      cat(sprintf("Release %s\n  against output root %s and stage root %s\n", o$release, output_root, stage_root))
      result <- compare_release(o$release, output_root, stage_root, tolerance, keys, allow)
    } else {
      if (length(o$positional) != 2) {
        compare_outputs_usage()
        return(2L)
      }
      cat(sprintf("Old %s\nNew %s\n", o$positional[1], o$positional[2]))
      prefix <- if (is.null(o$prefix)) "" else o$prefix
      result <- compare_output_trees(o$positional[1], o$positional[2], prefix, tolerance, keys, allow)
    }
    print_comparison(result)
    if (!is.null(o$out)) {
      write_comparison(result, o$out)
      cat(sprintf("Differences written to %s (%d rows)\n", o$out, nrow(result$diffs)))
    }
    if (result$unexpected > 0) 1L else 0L
  }, error = function(e) {
    message("compare_outputs: ", conditionMessage(e))
    2L
  })
}

if (sys.nframe() == 0L) quit(status = compare_outputs_main(), save = "no")
