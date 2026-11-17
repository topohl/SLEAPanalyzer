# ============================================================================
# 35_hidden_sensitivity.R
#
# Sensitivity of the Exp9 social preference results to the hidden-contact
# rule. Compares two SocP output trees scored from the same tracking and
# decisions, one with hidden time counted as contact (hidden_as: contact, the
# analysis) and one with it unobserved (hidden_as: unobserved: neither contact
# nor valid time), recording by recording. Descriptive only, no tests.
#
#   Rscript 35_hidden_sensitivity.R <contact_root> <unobserved_root> <out_dir> [--groups groups.csv]
#
#   *_root      <root>/<batch>/<phase>/combined_output.csv, as DLCA_SocP
#               writes it (output_dir of config/socp_all.yaml and
#               socp_hab_all.yaml)
#   out_dir     written: recordings.csv (batch, phase, Code, file, measure,
#               as_contact, as_unobserved, change = as_unobserved - as_contact)
#               and, with --groups, groups.csv (phase, measure, group, and n
#               and median under each variant) and group_differences.csv
#               (phase, measure, comparison "a - b" for each pair of groups in
#               the order they first appear, difference of the medians under
#               each variant)
#   --groups    a CSV with columns Code and group (any labels); animals with
#               no group are left out of the group tables and named
#
# Measures: contactLeft, contactRight, contactNovel, contactFamiliar,
# validTime, hiddenLeft and hiddenRight of combined_output.csv, and the
# preference index (contactNovel - contactFamiliar) / (contactNovel +
# contactFamiliar) as 09_assemble_all_batches.R computes it (missing without
# contact). Medians leave out missing values. Checks, each stopping the script
# with the recordings named: both trees hold the same recordings, each once;
# hidden time is the same in both (same tracking and decisions); and in every
# recording the contact tree's valid time exceeds the unobserved tree's by
# exactly the hidden time (the variants in this order). Runtime: seconds.
# ============================================================================

suppressMessages(library(data.table))
USAGE <- "usage: Rscript 35_hidden_sensitivity.R <contact_root> <unobserved_root> <out_dir> [--groups groups.csv]"
args <- commandArgs(TRUE)
GROUPS <- NULL
at <- which(args == "--groups")
if (length(at)) {
  if (length(at) > 1 || at == length(args)) stop(USAGE)
  GROUPS <- args[at + 1]; args <- args[-c(at, at + 1)]
}
if (length(args) != 3 || any(startsWith(args, "--"))) stop(USAGE)
ROOTS <- c(contact = args[1], unobserved = args[2]); OUT <- args[3]
for (r in ROOTS) if (!dir.exists(r)) stop("no SocP output tree at ", r)
if (normalizePath(ROOTS[[1]]) == normalizePath(ROOTS[[2]])) stop("the two trees are the same folder")
if (!is.null(GROUPS) && !file.exists(GROUPS)) stop("no groups file at ", GROUPS)
TOL <- 1e-6                                          # seconds; all times are whole frames

READ <- c("contactLeft", "contactRight", "contactNovel", "contactFamiliar", "validTime", "hiddenLeft", "hiddenRight")
MEASURES <- c(READ[1:4], "pref_index", READ[5:7])
read_tree <- function(root) {
  batches <- grep("^B[0-9]+$", basename(list.dirs(root, recursive = FALSE)), value = TRUE)
  files <- unlist(lapply(batches, function(b)
    file.path(root, b, basename(list.dirs(file.path(root, b), recursive = FALSE)), "combined_output.csv")))
  files <- files[file.exists(files)]
  if (!length(files)) stop("no SocP output (<batch>/<phase>/combined_output.csv) under ", root)
  rbindlist(lapply(files, function(f) {
    d <- fread(f, colClasses = list(character = c("file", "Code")))
    lacking <- setdiff(c("file", "Code", READ), names(d))
    if (length(lacking)) stop(f, " has no column ", paste(lacking, collapse = ", "))
    d <- d[, c("file", "Code", READ), with = FALSE]
    d[, pref_index := (contactNovel - contactFamiliar) / (contactNovel + contactFamiliar)]
    d[!is.finite(pref_index), pref_index := NA_real_]
    d[, `:=`(batch = basename(dirname(dirname(f))), phase = basename(dirname(f)))]
    setcolorder(d, c("batch", "phase"))
  }))
}
trees <- lapply(ROOTS, read_tree)
KEY <- c("batch", "phase", "Code", "file")
name_rows <- function(d) paste(paste(d$batch, d$phase, d$file, sep = "/"), collapse = ", ")
for (v in names(trees)) {
  twice <- trees[[v]][duplicated(trees[[v]][, ..KEY])]
  if (nrow(twice)) stop("recordings listed twice in the ", v, " tree: ", name_rows(twice))
}
only <- lapply(names(trees), function(v) fsetdiff(trees[[v]][, ..KEY], trees[[setdiff(names(trees), v)]][, ..KEY]))
if (any(vapply(only, nrow, 0L) > 0)) {
  stop("the trees hold different recordings; only in the contact tree: ", name_rows(only[[1]]),
       "; only in the unobserved tree: ", name_rows(only[[2]]))
}
both <- merge(trees$contact, trees$unobserved, by = KEY, suffixes = c(".contact", ".unobserved"))

differ <- both[abs(hiddenLeft.contact - hiddenLeft.unobserved) > TOL | abs(hiddenRight.contact - hiddenRight.unobserved) > TOL]
if (nrow(differ)) stop("hidden time differs between the trees (not scored from the same tracking and decisions): ", name_rows(differ))
if (!any(both$hiddenLeft.contact + both$hiddenRight.contact > 0)) {
  stop("no recording has hidden time: were the trees scored with hidden_contact: true?")
}
off <- both[abs(validTime.contact - validTime.unobserved - (hiddenLeft.contact + hiddenRight.contact)) > TOL]
if (nrow(off)) {
  stop("valid time under contact minus under unobserved is not the hidden time (trees swapped, or not scored ",
       "with hidden_as contact and unobserved?): ", name_rows(off))
}

long <- rbindlist(lapply(MEASURES, function(m) both[, .(batch, phase, Code, file, measure = m,
  as_contact = get(paste0(m, ".contact")), as_unobserved = get(paste0(m, ".unobserved")))]))
long[, change := as_unobserved - as_contact]
long <- long[order(as.integer(sub("^B", "", batch)), phase, Code, match(measure, MEASURES))]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fwrite(long, file.path(OUT, "recordings.csv"))

hidden <- both[hiddenLeft.contact + hiddenRight.contact > 0]
cat(sprintf("SocP recordings: %d; with hidden time: %d (%.1f s in all); the preference index by phase:\n",
            nrow(both), nrow(hidden), sum(hidden$hiddenLeft.contact + hidden$hiddenRight.contact)))
print(long[measure == "pref_index", .(recordings = .N, changed = sum(abs(change) > TOL, na.rm = TRUE),
                                      largest_change = if (all(is.na(change))) NA_real_ else max(abs(change), na.rm = TRUE)),
           by = phase], row.names = FALSE, class = FALSE)
written <- "recordings.csv"

if (!is.null(GROUPS)) {
  grp <- fread(GROUPS, colClasses = "character")
  if (!all(c("Code", "group") %in% names(grp))) stop(GROUPS, " needs the columns Code and group")
  grp <- unique(grp[!is.na(group) & nzchar(group), .(Code, group)])
  twice <- grp[duplicated(Code)]
  if (nrow(twice)) stop("more than one group for: ", paste(unique(twice$Code), collapse = ", "))
  again <- both[duplicated(both[, .(phase, Code)])]
  if (nrow(again)) stop("an animal recorded twice in a phase cannot be counted once: ", name_rows(again))
  ungrouped <- setdiff(both$Code, grp$Code)
  if (length(ungrouped)) cat("no group in", GROUPS, "(left out of the group tables):", paste(sort(ungrouped), collapse = ", "), "\n")
  group_order <- unique(grp$group)
  grouped <- grp[long, on = "Code", nomatch = NULL]
  med <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
  groups <- grouped[, .(n_as_contact = sum(!is.na(as_contact)), median_as_contact = med(as_contact),
                  n_as_unobserved = sum(!is.na(as_unobserved)), median_as_unobserved = med(as_unobserved)),
              by = .(phase, measure, group)]
  groups <- groups[order(phase, match(measure, MEASURES), match(group, group_order))]
  pairs <- if (length(group_order) > 1) as.data.table(t(utils::combn(group_order, 2))) else data.table(V1 = character(), V2 = character())
  pairs[, comparison := paste(V1, "-", V2)]
  differences <- rbindlist(lapply(seq_len(nrow(pairs)), function(k) {
    a <- groups[group == pairs$V1[k]]; b <- groups[group == pairs$V2[k]]
    ab <- merge(a, b, by = c("phase", "measure"), suffixes = c(".a", ".b"))
    ab[, .(phase, measure, comparison = pairs$comparison[k],
           as_contact = median_as_contact.a - median_as_contact.b,
           as_unobserved = median_as_unobserved.a - median_as_unobserved.b)]
  }))
  if (nrow(differences)) differences <- differences[order(phase, match(measure, MEASURES), match(comparison, pairs$comparison))]
  fwrite(groups, file.path(OUT, "groups.csv"))
  fwrite(differences, file.path(OUT, "group_differences.csv"))
  written <- c(written, "groups.csv", "group_differences.csv")
}
cat("wrote", paste(file.path(OUT, written), collapse = ", "), "\n")
