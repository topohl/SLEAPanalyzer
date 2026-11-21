#!/usr/bin/env Rscript
# Rescore one assay batch (social preference: one batch and phase) from a given
# code checkout, with a backup of the output it replaces and checks afterwards.
#
# Usage:
#   Rscript tools/rescore_batch.R --assay NOR --batch B1                          # plan only
#   Rscript tools/rescore_batch.R --assay NOR --batch B1 --tag borrowed_frames --allow allow.csv --write
#   Options: --phase (SocP: HAB, S1 or S2; NOR: NOV, the default, or HAB), --repo CHECKOUT (default:
#   the repository holding this tool), --config FILE (default: the checkout's configuration below),
#   --tag TAG (default: the checkout's short commit), --rscript PATH, --tolerance 1e-9.
#   Or source() this file and call rescore_plan(), print_rescore_plan() and run_rescore().
#
# Scripts (02_SLEAPanalzyer/) and configurations (validation/exp9_boris/config/) of the checkout:
#   SocP  DLCA_SocP v.1.1.0.R  socp_all.yaml (S1, S2), socp_hab_all.yaml (HAB)
#   EPM   DLCA_EPM v1.0.0.R    epm_all.yaml
#   NOR   DLCA_NOR v1.2.1.R    nor_all.yaml (NOV), nor_hab_all.yaml (HAB)
#   OFT   DLCA_OFT v1.2.0.R    oft_all.yaml
#
# Without --write it prints the plan and writes nothing. With --write:
#   1. It refuses when the checkout has uncommitted changes to tracked files or,
#      for the assays that write a run manifest (all but OFT), untracked files:
#      the manifest would then record dirty_working_tree: yes.
#   2. It copies the batch's live output folder (<output_dir>/<B>, SocP
#      <output_dir>/<B>/<phase>, OFT <behavior_root>/<B>/<output_subdir>) to
#      _superseded_<date>_<tag>/<label>_before/ next to the output tree (OFT: next
#      to the staging tree's OFT folder), as tools/retrack_recording.py does, with
#      an md5 manifest <label>_before_md5.csv, and checks the copy. It never
#      overwrites: a second run with the same date and tag takes <label>__2.
#   3. It writes <label>.yaml there: the assay configuration with batches (and
#      phases) set to this batch and every relative path made absolute against
#      the configuration's own folder, as the loader resolves it. The path keys
#      are the ones the checkout's loader resolves (rescore_path_keys()). It
#      checks that the copy reads as the original apart from those values. The
#      plan warns about a path that lies inside another git checkout.
#   4. It runs the checkout's assay script with SLEAP_ANALYZER_CONFIG pointing at
#      that copy (log <label>_run.log). OFT also rewrites its all-batches folder
#      OFT_SLEAP_output_v1.2.0 with this batch alone, so that folder is backed up
#      first and put back afterwards; the one-batch version is kept as
#      <label>_all_batches_written/. Its all_batches_*.csv include the rescored
#      batch only after a full OFT run.
#   5. It checks that the batch's run_manifest.yaml was written by this run and
#      names the checkout's HEAD, dirty_working_tree: no, this batch and phase,
#      the configuration copy and the decisions file the copy names (path and
#      md5). OFT writes no manifest: its checkout state is the one checked in
#      step 1, recorded in the log, and a failures.csv from this run fails it.
#      The files of the output folder the run did not rewrite are listed, for
#      information, here and in the log.
#   6. It compares the new output with the backup (tools/compare_outputs.R, file
#      paths as in a release: NOR/B1/combined_output.csv), with --allow as the
#      allowlist, and writes every difference to <label>_compare.csv.
#
# Exit status: 0 when every check passes and no difference is unexpected (plan
# only: when nothing stops the run), 1 otherwise, 2 usage or input error.

# The folder of this file: the source() call reading it, else the script Rscript runs.
# (When a script run by Rscript sources this file, --file names that script.)
rescore_tool_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(dirname(normalizePath(f, winslash = "/")))
  }
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  normalizePath("tools", winslash = "/", mustWork = FALSE)
}

RESCORE_TOOL_DIR <- rescore_tool_dir()

if (!exists("compare_output_trees", mode = "function")) {
  source(file.path(RESCORE_TOOL_DIR, "compare_outputs.R"), local = environment())
}

# The assay scripts and their configurations by phase ("" where an assay has none).
rescore_assays <- function() {
  list(
    SocP = list(script = "DLCA_SocP v.1.1.0.R",
                configs = c(HAB = "socp_hab_all.yaml", S1 = "socp_all.yaml", S2 = "socp_all.yaml")),
    EPM = list(script = "DLCA_EPM v1.0.0.R", configs = c("epm_all.yaml")),
    NOR = list(script = "DLCA_NOR v1.2.1.R", configs = c(NOV = "nor_all.yaml", HAB = "nor_hab_all.yaml")),
    OFT = list(script = "DLCA_OFT v1.2.0.R", configs = c("oft_all.yaml"))
  )
}

#' The configuration keys a checkout's loader resolves against the configuration's folder.
#'
#' Read from the checkout's code, so a path field added there is picked up: for EPM,
#' NOR and SocP the path_fields load_assay_config() (core/assay_config.R) hands to
#' load_analysis_config(), for OFT the path_fields of the apply_config_overlay() call
#' in its script.
rescore_path_keys <- function(assay, repo) {
  code <- file.path(repo, "02_SLEAPanalzyer")
  if (assay == "OFT") {
    script <- file.path(code, rescore_assays()$OFT$script)
    if (!file.exists(script)) stop("the checkout has no ", file.path("02_SLEAPanalzyer", basename(script)))
    found <- list()
    walk <- function(e) {
      if (!is.call(e)) return(invisible())
      if (identical(e[[1]], as.name("apply_config_overlay")) && "path_fields" %in% names(e)) {
        found[[length(found) + 1]] <<- eval(e[["path_fields"]], baseenv())
      }
      for (i in seq_along(e)) if (is.call(e[[i]])) walk(e[[i]])
    }
    for (e in parse(script, keep.source = FALSE)) walk(e)
    keys <- unique(found)
    if (length(keys) != 1) {
      stop("cannot read the path fields of ", basename(script), ": ", length(keys),
           " different apply_config_overlay(path_fields = ...) calls")
    }
    return(as.character(keys[[1]]))
  }
  files <- file.path(code, "core", c("config.R", "assay_config.R"))
  if (!all(file.exists(files))) stop("the checkout has no ", paste(files[!file.exists(files)], collapse = ", "))
  env <- new.env(parent = baseenv())
  for (f in files) sys.source(f, envir = env)
  # The loader's own path_fields, without reading any file.
  env$load_analysis_config <- function(path, defaults = list(), schema = NULL, path_fields = character()) path_fields
  keys <- tryCatch(env$load_assay_config("(none)", assay), error = function(e) {
    stop("cannot read the path fields of ", files[2], ": ", conditionMessage(e), call. = FALSE)
  })
  if (!is.character(keys)) stop("cannot read the path fields of ", files[2])
  keys
}

# TRUE where two paths name the same file or folder.
same_path <- function(a, b) {
  norm <- function(x) normalizePath(path.expand(x), winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows") tolower(norm(a)) == tolower(norm(b)) else norm(a) == norm(b)
}

# The top folder of the git checkout holding a path (or its nearest existing parent), else NULL.
checkout_of <- function(path) {
  d <- path.expand(path)
  while (!dir.exists(d)) {
    parent <- dirname(d)
    if (identical(parent, d)) return(NULL)
    d <- parent
  }
  top <- git_lines(d, "rev-parse", "--show-toplevel")
  if (length(top) == 0 || !nzchar(top[1])) NULL else top[1]
}

# The loader's test for an absolute path (core/config.R).
is_absolute_config_path <- function(x) grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", x)

default_rscript <- function() {
  file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
}

# ---------------------------------------------------------------------------
# The configuration copy
# ---------------------------------------------------------------------------

# The value after "key:" on a top-level line: list(value, quote, rest), where
# rest is what follows the value (a trailing comment, with its spacing).
split_yaml_value <- function(text) {
  s <- sub("^[ \t]+", "", text)
  if (startsWith(s, "\"")) {
    m <- regmatches(s, regexpr("^\"([^\"\\\\]|\\\\.)*\"", s, perl = TRUE))
    if (length(m) == 0) stop("unclosed double quote in: ", text)
    value <- gsub("\\\\(.)", "\\1", substr(m, 2, nchar(m) - 1), perl = TRUE)
    return(list(value = value, quote = "\"", rest = substring(s, nchar(m) + 1)))
  }
  if (startsWith(s, "'")) {
    m <- regmatches(s, regexpr("^'([^']|'')*'", s))
    if (length(m) == 0) stop("unclosed single quote in: ", text)
    value <- gsub("''", "'", substr(m, 2, nchar(m) - 1), fixed = TRUE)
    return(list(value = value, quote = "'", rest = substring(s, nchar(m) + 1)))
  }
  hash <- regexpr("(^|[ \t])#", s)
  value <- sub("[ \t]+$", "", if (hash > 0) substr(s, 1, hash - 1) else s)
  list(value = value, quote = "", rest = substring(s, nchar(value) + 1))
}

# A value written back in its original quoting, quoted where YAML needs it.
format_yaml_value <- function(value, quote) {
  if (quote == "'") return(paste0("'", gsub("'", "''", value, fixed = TRUE), "'"))
  plain_ok <- quote == "" && !grepl("^[-?:,\\[\\]{}#&*!|>'\"%@`\\s]|: | #|:$|\\s$", value, perl = TRUE)
  if (plain_ok) return(value)
  paste0("\"", gsub("([\"\\\\])", "\\\\\\1", value, perl = TRUE), "\"")
}

# Replace a top-level key with one line (its continuation lines go too), or add it at the end.
set_yaml_line <- function(lines, key, line) {
  at <- grep(paste0("^", key, ":"), lines)
  if (length(at) > 1) stop("the configuration sets ", key, " more than once")
  if (length(at) == 0) return(c(lines, line))
  end <- at
  while (end < length(lines) && grepl("^([ \t]+[^ \t]|-)", lines[end + 1])) end <- end + 1
  c(lines[seq_len(at - 1)], line, lines[-seq_len(end)])
}

#' The configuration text with batches (and phases) set and relative paths made absolute.
#'
#' @param lines the configuration file's lines
#' @param set named list of keys to set to one value each, as a one-item list (batches = "B1")
#' @param path_keys keys holding paths
#' @param config_dir the folder relative paths are resolved against
#' @return list(lines, changes): changes is a data frame of key, old, new
rewrite_config_lines <- function(lines, set, path_keys, config_dir) {
  changes <- data.frame(key = character(), old = character(), new = character(), stringsAsFactors = FALSE)
  for (key in names(set)) {
    at <- grep(paste0("^", key, ":"), lines)
    old <- if (length(at)) trimws(sub(paste0("^", key, ":"), "", lines[at[1]])) else "(not set)"
    new <- paste0("[", paste(set[[key]], collapse = ", "), "]")
    lines <- set_yaml_line(lines, key, paste0(key, ": ", new))
    changes[nrow(changes) + 1, ] <- list(key, old, new)
  }
  for (key in path_keys) {
    at <- grep(paste0("^", key, ":"), lines)
    if (length(at) > 1) stop("the configuration sets ", key, " more than once")
    if (length(at) == 0) next
    v <- split_yaml_value(sub(paste0("^", key, ":"), "", lines[at]))
    if (!nzchar(v$value) || (v$quote == "" && tolower(v$value) %in% c("~", "null"))) next
    if (is_absolute_config_path(v$value)) next
    new <- normalizePath(file.path(config_dir, v$value), winslash = "/", mustWork = FALSE)
    lines[at] <- paste0(key, ": ", format_yaml_value(new, v$quote), v$rest)
    changes[nrow(changes) + 1, ] <- list(key, v$value, new)
  }
  list(lines = lines, changes = changes)
}

#' Problems with a configuration copy: it must read as the original except for the set keys
#' and the relative paths, which must read as the loader would resolve them.
check_config_copy <- function(original, copy, set, path_keys, config_dir) {
  a <- yaml::yaml.load(paste(original, collapse = "\n"))
  b <- yaml::yaml.load(paste(copy, collapse = "\n"))
  problems <- character()
  for (k in union(names(a), names(b))) {
    if (k %in% names(set)) {
      if (!identical(as.character(unlist(b[[k]])), as.character(set[[k]]))) {
        problems <- c(problems, sprintf("%s reads %s, not %s", k, paste(unlist(b[[k]]), collapse = ","),
                                        paste(set[[k]], collapse = ",")))
      }
      next
    }
    v <- a[[k]]
    if (k %in% path_keys && is.character(v) && length(v) == 1 && nzchar(v) && !is_absolute_config_path(v)) {
      want <- normalizePath(file.path(config_dir, v), winslash = "/", mustWork = FALSE)
      if (!identical(b[[k]], want)) problems <- c(problems, sprintf("%s reads %s, not %s", k, format(b[[k]]), want))
      next
    }
    if (!identical(v, b[[k]])) problems <- c(problems, sprintf("%s changed", k))
  }
  problems
}

# ---------------------------------------------------------------------------
# The plan
# ---------------------------------------------------------------------------

git_lines <- function(repo, ...) {
  out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(repo), ...), stdout = TRUE, stderr = FALSE)),
                  error = function(e) NULL)
  status <- attr(out, "status")
  if (is.null(out) && is.null(status)) return(NULL)
  if (!is.null(status) && status != 0) return(NULL)
  out
}

# The first of label, label__2, label__3, ... with no backup, config copy or md5 manifest yet.
free_run_id <- function(root, label) {
  n <- 1
  id <- label
  while (any(file.exists(file.path(root, paste0(id, c("_before", "_before_md5.csv", ".yaml", "_run.log")))))) {
    n <- n + 1
    id <- paste0(label, "__", n)
  }
  id
}

#' Work out everything a rescoring of one batch would do, without writing anything.
#'
#' @return a list; its element `problems` lists what stops --write
rescore_plan <- function(assay, batch, phase = NULL, repo = NULL, config = NULL, tag = NULL, allow = NULL,
                         rscript = default_rscript(), date = Sys.Date()) {
  assays <- rescore_assays()
  if (!assay %in% names(assays)) stop("--assay must be one of ", paste(names(assays), collapse = ", "))
  if (!grepl("^[A-Za-z0-9]+$", batch)) stop("--batch must be letters and digits, such as B1")
  spec <- assays[[assay]]
  phases <- names(spec$configs)
  if (is.null(phases)) {
    if (!is.null(phase)) stop(assay, " has no phases; leave out --phase")
    config_name <- spec$configs[[1]]
  } else {
    if (is.null(phase)) {
      if (assay == "NOR") phase <- "NOV" else stop(assay, " needs --phase (", paste(phases, collapse = ", "), ")")
    }
    if (!phase %in% phases) stop("--phase for ", assay, " must be one of ", paste(phases, collapse = ", "))
    config_name <- spec$configs[[phase]]
  }

  if (is.null(repo)) repo <- dirname(RESCORE_TOOL_DIR)
  if (!dir.exists(repo)) stop("checkout not found: ", repo)
  repo <- normalizePath(repo, winslash = "/")
  script <- file.path(repo, "02_SLEAPanalzyer", spec$script)
  if (!file.exists(script)) stop("the checkout has no ", file.path("02_SLEAPanalzyer", spec$script))
  if (is.null(config)) config <- file.path(repo, "validation", "exp9_boris", "config", config_name)
  if (!file.exists(config)) stop("configuration not found: ", config)
  config <- normalizePath(config, winslash = "/")
  config_dir <- dirname(config)
  original <- readLines(config, warn = FALSE, encoding = "UTF-8")
  parsed <- yaml::yaml.load(paste(original, collapse = "\n"))

  problems <- character()
  head <- git_lines(repo, "rev-parse", "HEAD")
  if (is.null(head) || length(head) == 0) {
    problems <- c(problems, paste(repo, "is not a git checkout"))
    head <- NA_character_
  }
  status <- git_lines(repo, "status", "--porcelain")
  if (is.null(status)) status <- character()
  tracked <- status[!startsWith(status, "??")]
  untracked <- status[startsWith(status, "??")]
  if (length(tracked) > 0) {
    problems <- c(problems, sprintf("the checkout has uncommitted changes to tracked files: %s",
                                    paste(utils::head(trimws(tracked), 5), collapse = "; ")))
  }
  writes_manifest <- assay != "OFT"
  if (writes_manifest && length(untracked) > 0) {
    problems <- c(problems, sprintf(paste(
      "the checkout has untracked files, so the run manifest would record dirty_working_tree: yes;",
      "commit, remove or ignore them: %s"), paste(utils::head(sub("^\\?\\? ", "", untracked), 5), collapse = "; ")))
  }
  if (writes_manifest && isFALSE(parsed$write_manifest)) {
    problems <- c(problems, paste(basename(config), "sets write_manifest: false; the manifest check needs it"))
  }
  if (!is.null(parsed$batches) && !batch %in% parsed$batches) {
    problems <- c(problems, sprintf("%s is not among the batches of %s (%s)", batch, basename(config),
                                    paste(parsed$batches, collapse = ", ")))
  }
  if (assay == "SocP" && !is.null(parsed$phases) && !phase %in% parsed$phases) {
    problems <- c(problems, sprintf("%s is not among the phases of %s (%s)", phase, basename(config),
                                    paste(parsed$phases, collapse = ", ")))
  }

  resolve <- function(v) {
    if (is.null(v) || !is.character(v) || !nzchar(v)) return(NULL)
    if (is_absolute_config_path(v)) v else normalizePath(file.path(config_dir, v), winslash = "/", mustWork = FALSE)
  }
  socp_phase <- if (assay == "SocP") phase else NULL
  global_dir <- NULL
  if (assay == "OFT") {
    root <- resolve(parsed$behavior_root)
    if (is.null(root)) stop(basename(config), " sets no behavior_root")
    input_sub <- if (is.null(parsed$input_subdir)) file.path("OFT", "formatted") else parsed$input_subdir
    output_sub <- if (is.null(parsed$output_subdir)) file.path("OFT", "output_v1.2.0") else parsed$output_subdir
    input_dir <- normalizePath(file.path(root, batch, input_sub), winslash = "/", mustWork = FALSE)
    out_dir <- file.path(root, batch, output_sub)
    global_dir <- file.path(root, "OFT_SLEAP_output_v1.2.0")
    backup_root <- file.path(dirname(root), paste0("_superseded_", date, "_"))
    prefix <- paste("OFT", batch, sep = "/")
  } else {
    input_root <- resolve(parsed$input_dir)
    output_root <- resolve(parsed$output_dir)
    if (is.null(input_root) || is.null(output_root)) stop(basename(config), " sets no input_dir or output_dir")
    input_dir <- do.call(file.path, as.list(c(input_root, batch, socp_phase)))
    out_dir <- do.call(file.path, as.list(c(output_root, batch, socp_phase)))
    backup_root <- file.path(dirname(dirname(output_root)), paste0("_superseded_", date, "_"))
    prefix <- paste(c(if (assay == "NOR" && phase == "HAB") "NOR_HAB" else assay, batch, socp_phase), collapse = "/")
  }
  if (is.null(tag)) tag <- if (is.na(head)) "rescore" else substr(head, 1, 7)
  if (!grepl("^[A-Za-z0-9._-]+$", tag)) stop("--tag must be letters, digits, dot, underscore or hyphen")
  backup_root <- paste0(backup_root, tag)
  label <- gsub("/", "_", prefix, fixed = TRUE)
  id <- free_run_id(backup_root, label)
  at <- function(suffix) file.path(backup_root, paste0(id, suffix))

  set <- list(batches = batch)
  if (assay == "SocP") set$phases <- phase
  path_keys <- rescore_path_keys(assay, repo)
  copy <- rewrite_config_lines(original, set, path_keys, config_dir)
  copy_problems <- check_config_copy(original, copy$lines, set, path_keys, config_dir)
  if (length(copy_problems) > 0) {
    problems <- c(problems, paste("the configuration copy does not read as intended:", paste(copy_problems, collapse = "; ")))
  }
  # A path inside another checkout reads that checkout's files, not this one's.
  copied <- yaml::yaml.load(paste(copy$lines, collapse = "\n"))
  warnings <- character()
  for (k in path_keys) {
    v <- copied[[k]]
    if (!is.character(v) || length(v) != 1 || !nzchar(v)) next
    top <- checkout_of(v)
    if (!is.null(top) && !same_path(top, repo)) {
      warnings <- c(warnings, sprintf("%s lies inside another checkout (%s): %s", k, top, v))
    }
  }
  copied_file <- function(key) {
    v <- copied[[key]]
    if (is.character(v) && length(v) == 1 && nzchar(v)) v else NULL
  }
  decisions <- copied_file("qc_decisions_file")
  hand_events <- copied_file("hand_events_file")
  hidden_review <- copied_file("hidden_review_file")

  input_files <- if (dir.exists(input_dir)) length(list.files(input_dir, pattern = "\\.csv$")) else NA_integer_
  if (is.na(input_files) || input_files == 0) problems <- c(problems, paste("no tracking files in", input_dir))
  out_files <- if (dir.exists(out_dir)) list.files(out_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE) else character()
  global_files <- if (!is.null(global_dir) && dir.exists(global_dir)) {
    list.files(global_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  } else character()
  # R cannot open paths of 260 characters or more on this system.
  longest <- max(0, nchar(file.path(at("_before"), out_files)), nchar(file.path(at("_all_batches_before"), global_files)),
                 nchar(file.path(at("_all_batches_written"), global_files)))
  if (longest >= 260) problems <- c(problems, sprintf("a backup path would be %d characters long (limit 259)", longest))

  if (!is.null(allow)) {
    allow <- normalizePath(allow, winslash = "/", mustWork = FALSE)
    checked <- tryCatch(read_allowlist(allow), error = function(e) conditionMessage(e))
    if (is.character(checked)) problems <- c(problems, checked)
  }
  if (!file.exists(rscript)) problems <- c(problems, paste("Rscript not found:", rscript))

  list(
    assay = assay, batch = batch, phase = phase, socp_phase = socp_phase, label = label, prefix = prefix,
    repo = repo, head = head, tracked = tracked, untracked = untracked, script = script, rscript = rscript,
    config = config, config_lines = original, copy_lines = copy$lines, changes = copy$changes,
    path_keys = path_keys, decisions = decisions, hand_events = hand_events, hidden_review = hidden_review,
    input_dir = input_dir, input_files = input_files, out_dir = out_dir, out_files = length(out_files),
    global_dir = global_dir, writes_manifest = writes_manifest,
    backup_root = backup_root, id = id, backup_dir = at("_before"), backup_md5 = at("_before_md5.csv"),
    global_backup = at("_all_batches_before"), global_md5 = at("_all_batches_before_md5.csv"),
    global_written = at("_all_batches_written"), config_copy = at(".yaml"), log = at("_run.log"),
    compare_out = at("_compare.csv"), allow = allow, warnings = warnings, problems = problems
  )
}

print_rescore_plan <- function(plan) {
  state <- if (length(plan$tracked)) sprintf("%d uncommitted changes", length(plan$tracked)) else "clean"
  if (length(plan$untracked)) state <- sprintf("%s, %d untracked files", state, length(plan$untracked))
  cat(sprintf("Rescore %s\n", paste(c(plan$assay, plan$batch, plan$phase), collapse = " ")))
  cat(sprintf("  checkout  %s at %s (%s)\n", plan$repo, plan$head, state))
  cat(sprintf("  script    %s\n", basename(plan$script)))
  cat(sprintf("  config    %s\n", plan$config))
  cat(sprintf("            copy %s with:\n", plan$config_copy))
  for (i in seq_len(nrow(plan$changes))) {
    cat(sprintf("              %s: %s -> %s\n", plan$changes$key[i], plan$changes$old[i], plan$changes$new[i]))
  }
  cat(sprintf("  input     %s (%s tracking files)\n", plan$input_dir, plan$input_files))
  cat(sprintf("  output    %s (%d files%s)\n", plan$out_dir, plan$out_files,
              if (plan$out_files == 0) "; nothing to back up" else ""))
  cat(sprintf("  backup    %s, md5 manifest %s\n", plan$backup_dir, basename(plan$backup_md5)))
  if (!is.null(plan$global_dir)) {
    cat(sprintf("  all batches  %s is rewritten with %s alone: backed up to %s and put back afterwards\n",
                plan$global_dir, plan$batch, basename(plan$global_backup)))
  }
  cat(sprintf("  log       %s\n", plan$log))
  if (plan$writes_manifest) {
    cat(sprintf("  check     run_manifest.yaml names %s, dirty_working_tree: no, batch %s%s,\n", plan$head,
                plan$batch, if (!is.null(plan$socp_phase)) paste(", phase", plan$socp_phase) else ""))
    cat(sprintf("            config %s and %s\n", basename(plan$config_copy),
                if (is.null(plan$decisions)) "no qc_decisions_file" else
                  sprintf("qc_decisions_file %s (path and md5)", plan$decisions)))
    if (!is.null(plan$hand_events)) cat(sprintf("            and hand_events_file %s (path and md5)\n", plan$hand_events))
    if (!is.null(plan$hidden_review)) cat(sprintf("            and hidden_review_file %s (path and md5)\n", plan$hidden_review))
  } else {
    cat("  check     OFT writes no run manifest: the checkout state above goes into the log; a failures.csv fails\n")
  }
  cat(sprintf("  compare   new output against the backup%s; differences to %s\n",
              if (is.null(plan$allow)) "" else paste(" with allowlist", plan$allow), basename(plan$compare_out)))
  if (length(plan$warnings)) {
    cat("Warning:\n")
    cat(paste0("  - ", plan$warnings, "\n"), sep = "")
  }
  if (length(plan$problems)) {
    cat("Will not run:\n")
    cat(paste0("  - ", plan$problems, "\n"), sep = "")
  }
  invisible(plan)
}

# ---------------------------------------------------------------------------
# Doing it
# ---------------------------------------------------------------------------

# Copy every file under src to dst (which must not exist yet); returns the relative paths.
copy_tree <- function(src, dst) {
  if (file.exists(dst)) stop("will not overwrite ", dst)
  files <- list.files(src, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  dir.create(dst, recursive = TRUE)
  for (d in unique(dirname(files[dirname(files) != "."]))) dir.create(file.path(dst, d), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(file.path(src, files), file.path(dst, files), overwrite = FALSE, copy.date = TRUE)
  if (!all(ok)) stop("could not copy to ", dst, ": ", paste(utils::head(files[!ok], 5), collapse = ", "))
  files
}

#' Copy a folder, check every copy against its source by md5 and write the md5 manifest.
backup_folder <- function(src, dst, manifest) {
  if (file.exists(manifest)) stop("will not overwrite ", manifest)
  files <- if (dir.exists(src)) copy_tree(src, dst) else {
    if (file.exists(dst)) stop("will not overwrite ", dst)
    dir.create(dst, recursive = TRUE)
    character()
  }
  md5_src <- unname(tools::md5sum(file.path(src, files)))
  md5_dst <- unname(tools::md5sum(file.path(dst, files)))
  if (any(md5_src != md5_dst)) stop("the backup differs from its source: ", paste(files[md5_src != md5_dst], collapse = ", "))
  utils::write.csv(data.frame(file = files, bytes = file.size(file.path(dst, files)), md5 = md5_dst),
                   manifest, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(files)
}

#' Checks of a run manifest: written since `since`, the checkout's HEAD, a clean tree, this batch and
#' phase, the configuration copy (`config_copy`) and, where the manifest records its decisions, the
#' qc_decisions_file the copy names (`decisions`, NULL for none) by path and md5.
#'
#' @return data frame of check, ok, detail
check_run_manifest <- function(path, head, batch, phase = NULL, since = NULL, config_copy = NULL, decisions = NULL,
                               hand_events = NULL, hidden_review = NULL) {
  out <- data.frame(check = character(), ok = logical(), detail = character(), stringsAsFactors = FALSE)
  add <- function(check, ok, detail) out[nrow(out) + 1, ] <<- list(check, isTRUE(ok), paste(detail, collapse = " "))
  if (!file.exists(path)) {
    add("run manifest written", FALSE, paste("no", path))
    return(out)
  }
  m <- yaml::read_yaml(path)
  if (!is.null(since)) add("written by this run", file.mtime(path) >= since - 2, format(file.mtime(path)))
  commit <- if (is.null(m$code$commit)) "none" else as.character(m$code$commit)
  add("commit is the checkout's HEAD", identical(commit, head), commit)
  dirty <- m$code$dirty_working_tree
  add("dirty_working_tree: no", isFALSE(dirty), if (is.null(dirty)) "not recorded" else if (isTRUE(dirty)) "yes" else "no")
  got <- if (is.null(m$extra$batch)) "none" else as.character(m$extra$batch)
  add(paste("batch", batch), identical(got, batch), got)
  if (!is.null(phase)) {
    got <- if (is.null(m$extra$phase)) "none" else as.character(m$extra$phase)
    add(paste("phase", phase), identical(got, phase), got)
  }
  if (!is.null(config_copy)) {
    got <- m$config$.config_path
    read <- is.character(got) && length(got) == 1
    add("config is the copy", read && same_path(got, config_copy), if (read) got else "not recorded")
  }
  # A file the run read, as the manifest records it under extra$<name>, against the one the copy names.
  recorded_file <- function(name, key, want) {
    if (!name %in% names(m$extra)) return(invisible())
    rec <- m$extra[[name]]
    if (is.null(want)) {
      add(sprintf("no %s, as in the copy", key), length(rec) == 0,
          if (length(rec)) vapply(rec, function(r) format(r$path), "") else "none")
      return(invisible())
    }
    r <- if (length(rec) == 1) rec[[1]] else list()
    ok <- isTRUE(r$exists) && is.character(r$path) && same_path(r$path, want) &&
      identical(as.character(r$md5), unname(tools::md5sum(want)))
    add(sprintf("%s is the copy's (path, md5)", key), ok,
        if (length(rec) == 1) sprintf("%s, md5 %s", format(r$path), format(r$md5)) else sprintf("%d files", length(rec)))
  }
  recorded_file("qc_decisions", "qc_decisions_file", decisions)
  recorded_file("hand_events", "hand_events_file", hand_events)
  recorded_file("hidden_review", "hidden_review_file", hidden_review)
  out
}

# Run the assay script with SLEAP_ANALYZER_CONFIG set; stdout and stderr go to the log. Returns the exit status.
run_assay_script <- function(plan) {
  previous <- Sys.getenv("SLEAP_ANALYZER_CONFIG", unset = NA_character_)
  Sys.setenv(SLEAP_ANALYZER_CONFIG = plan$config_copy)
  old_wd <- setwd(plan$repo)
  on.exit({
    setwd(old_wd)
    if (is.na(previous)) Sys.unsetenv("SLEAP_ANALYZER_CONFIG") else Sys.setenv(SLEAP_ANALYZER_CONFIG = previous)
  }, add = TRUE)
  status <- suppressWarnings(system2(plan$rscript, c("--vanilla", shQuote(plan$script)), stdout = plan$log, stderr = plan$log))
  as.integer(status)
}

# OFT: keep the all-batches folder the one-batch run wrote, and put the previous one back.
restore_all_batches <- function(plan) {
  g <- plan$global_dir
  if (dir.exists(g)) {
    if (!suppressWarnings(file.rename(g, plan$global_written))) {
      tryCatch({
        copy_tree(g, plan$global_written)
        unlink(g, recursive = TRUE)
      }, error = function(e) NULL)
    }
    # A file held open (all_batches_summary.csv in a spreadsheet program) blocks both.
    if (dir.exists(g)) {
      stop(sprintf("%s could not be moved aside (is a file in it open?). Close the file, then by hand %s.", g,
                   if (dir.exists(plan$global_backup)) {
                     sprintf("replace that folder with the copy made before the run, %s", plan$global_backup)
                   } else {
                     "delete that folder; there was none before the run"
                   }), call. = FALSE)
    }
    cat(sprintf("  kept      the %s-only %s as %s\n", plan$batch, basename(g), plan$global_written))
  }
  if (dir.exists(plan$global_backup)) {
    copy_tree(plan$global_backup, g)
    cat(sprintf("  put back  %s\n", g))
  }
}

# Files under dir older than the run's start, so not rewritten by it (with the
# manifest check's 2 s margin for coarse file times).
not_rewritten <- function(dir, start) {
  files <- if (dir.exists(dir)) list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE) else character()
  files[file.mtime(file.path(dir, files)) < start - 2]
}

# The run log: a header, the assay script's output (already in the log), the checks and,
# after a successful run, the files of the output folder it did not rewrite (stale).
write_rescore_log <- function(plan, start, status, checks, stale) {
  output <- if (file.exists(plan$log)) readLines(plan$log, warn = FALSE) else character()
  writeLines(c(
    sprintf("rescore_batch.R %s", format(start, "%Y-%m-%d %H:%M:%S")),
    sprintf("checkout: %s", plan$repo),
    sprintf("commit: %s", plan$head),
    sprintf("uncommitted changes to tracked files: %s", if (length(plan$tracked)) length(plan$tracked) else "none"),
    sprintf("untracked files: %s", if (length(plan$untracked)) length(plan$untracked) else "none"),
    sprintf("script: %s", plan$script),
    sprintf("config: %s (copy of %s)", plan$config_copy, plan$config),
    if (length(plan$warnings)) paste("warning:", plan$warnings),
    sprintf("exit status: %d", status),
    "----",
    output,
    if (nrow(checks)) c("----", sprintf("%s %s (%s)", ifelse(checks$ok, "ok", "FAILED"), checks$check, checks$detail)),
    if (!is.null(stale)) c(sprintf("not rewritten by this run: %s", if (length(stale)) length(stale) else "none"),
                           if (length(stale)) paste0("  ", stale))
  ), plan$log)
}

#' Back up, run, check and compare. Returns 0 when every check passes and no difference is unexpected, else 1.
run_rescore <- function(plan, tolerance = 1e-9) {
  if (length(plan$problems)) stop("will not run: ", paste(plan$problems, collapse = "; "))
  dir.create(plan$backup_root, recursive = TRUE, showWarnings = FALSE)
  backup_folder(plan$out_dir, plan$backup_dir, plan$backup_md5)
  cat(sprintf("  backed up %s -> %s (%d files, md5 checked)\n", plan$out_dir, plan$backup_dir, plan$out_files))
  if (!is.null(plan$global_dir) && dir.exists(plan$global_dir)) {
    backup_folder(plan$global_dir, plan$global_backup, plan$global_md5)
    cat(sprintf("  backed up %s -> %s\n", plan$global_dir, plan$global_backup))
  }
  writeLines(enc2utf8(plan$copy_lines), plan$config_copy, useBytes = TRUE)

  cat(sprintf("  running   %s for %s ...\n", basename(plan$script), gsub("/", " ", sub("^[^/]+/", "", plan$prefix))))
  start <- Sys.time()
  status <- NA_integer_
  checks <- data.frame(check = character(), ok = logical(), detail = character(), stringsAsFactors = FALSE)
  stale <- NULL
  # The log gets its header however this ends, also when putting OFT's folder back stops it.
  on.exit(write_rescore_log(plan, start, status, checks, stale), add = TRUE)
  status <- tryCatch(run_assay_script(plan), error = function(e) {
    cat(sprintf("  could not start %s: %s\n", plan$rscript, conditionMessage(e)))
    -1L
  })
  if (!is.null(plan$global_dir)) {
    failures <- file.path(plan$global_dir, "failures.csv")
    if (file.exists(failures) && file.mtime(failures) >= start - 2) {
      f <- utils::read.csv(failures, colClasses = "character")
      mine <- f$file[f$Batch == plan$batch]
      checks[nrow(checks) + 1, ] <- list("no failed files", length(mine) == 0, paste(mine, collapse = ", "))
    }
    restore_all_batches(plan)
  }
  if (status != 0) {
    output <- if (file.exists(plan$log)) readLines(plan$log, warn = FALSE) else character()
    cat(sprintf("  FAILED    %s exited with status %d; last lines of %s:\n", basename(plan$script), status, plan$log))
    cat(paste0("            ", utils::tail(output, 15), "\n"), sep = "")
    cat(sprintf("  The output before the run is in %s.\n", plan$backup_dir))
    return(1L)
  }

  if (plan$writes_manifest) {
    checks <- rbind(checks, check_run_manifest(file.path(plan$out_dir, "run_manifest.yaml"), plan$head, plan$batch,
                                               plan$socp_phase, since = start, config_copy = plan$config_copy,
                                               decisions = plan$decisions, hand_events = plan$hand_events,
                                               hidden_review = plan$hidden_review))
  }
  for (i in seq_len(nrow(checks))) {
    cat(sprintf("  %-9s %s (%s)\n", if (checks$ok[i]) "ok" else "FAILED", checks$check[i], checks$detail[i]))
  }
  # Informational: such a file is in the backup unchanged, so it compares as identical.
  stale <- not_rewritten(plan$out_dir, start)
  if (length(stale)) {
    cat(sprintf("  note      %d %s of the output folder %s not rewritten by this run (older than its start):\n",
                length(stale), if (length(stale) == 1) "file" else "files", if (length(stale) == 1) "was" else "were"))
    cat(paste0("            ", utils::head(stale, 10), "\n"), sep = "")
    if (length(stale) > 10) cat(sprintf("            ... %d more in %s\n", length(stale) - 10, basename(plan$log)))
  }

  result <- compare_output_trees(plan$backup_dir, plan$out_dir, prefix = plan$prefix, tolerance = tolerance,
                                 allow = plan$allow)
  write_comparison(result, plan$compare_out)
  print_comparison(result)
  cat(sprintf("Differences written to %s (%d rows)\n", plan$compare_out, nrow(result$diffs)))
  if (!is.null(plan$global_dir)) {
    cat(sprintf("Note: %s/all_batches_*.csv do not include the rescored %s until a full OFT run.\n",
                plan$global_dir, plan$batch))
  }
  passed <- all(checks$ok) && result$unexpected == 0
  cat(if (passed) "Rescore done; every check passed.\n" else "Rescore done, but a check FAILED (see above).\n")
  if (passed) 0L else 1L
}

# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

rescore_batch_usage <- function() {
  cat(paste(
    "Usage: Rscript tools/rescore_batch.R --assay SocP|EPM|NOR|OFT --batch B1 [--phase PHASE]",
    "         [--repo CHECKOUT] [--config FILE] [--tag TAG] [--allow allow.csv] [--rscript PATH]",
    "         [--tolerance 1e-9] [--write]",
    "Without --write it prints the plan and writes nothing.",
    sep = "\n"), "\n")
}

rescore_batch_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  tryCatch({
    o <- parse_tool_args(argv, values = c("assay", "batch", "phase", "repo", "config", "tag", "allow", "rscript",
                                          "tolerance"), flags = c("write", "help"))
    if (isTRUE(o$help)) {
      rescore_batch_usage()
      return(0L)
    }
    if (is.null(o$assay) || is.null(o$batch) || length(o$positional) > 0) {
      rescore_batch_usage()
      return(2L)
    }
    tolerance <- if (is.null(o$tolerance)) 1e-9 else as.numeric(o$tolerance)
    if (!is.finite(tolerance) || tolerance < 0) stop("--tolerance must be a number >= 0")
    plan <- rescore_plan(o$assay, o$batch, phase = o$phase, repo = o$repo, config = o$config, tag = o$tag,
                         allow = o$allow, rscript = if (is.null(o$rscript)) default_rscript() else o$rscript)
    print_rescore_plan(plan)
    if (!isTRUE(o$write)) {
      cat("Plan only: nothing was written. Add --write to rescore.\n")
      return(if (length(plan$problems)) 1L else 0L)
    }
    if (length(plan$problems)) {
      cat("Nothing was written.\n")
      return(1L)
    }
    run_rescore(plan, tolerance)
  }, error = function(e) {
    message("rescore_batch: ", conditionMessage(e))
    2L
  })
}

if (sys.nframe() == 0L) quit(status = rescore_batch_main(), save = "no")
