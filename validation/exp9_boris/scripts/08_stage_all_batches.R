# ============================================================================
# 08_stage_all_batches.R
#
# Stages the SLEAP tracking exports for B1-B6 x {EPM, NOR, SocP, OFT} into a
# canonical layout the pipeline can consume without per-batch surprises.
# Reads Raw Data/ only; nothing there is modified.
#
# Staging lives OUTSIDE the iCloud project folder: it is ~6 GB of regenerable
# intermediate data and does not belong in a synced directory. Only pipeline
# OUTPUT goes back into the project.
#
# What this fixes, per 07_audit_batches.R:
#
#   1. ANIMAL CODE POSITION. EPM B2/B3/B4, OFT B1/B2 and the SocP B5 root
#      copies put the code somewhere other than the first four characters, so
#      the pipeline's str_extract(name, "^[A-Za-z0-9]{4}") returns NA and the
#      animal silently loses its metadata. Every staged file is renamed
#      <CODE>... so the rule holds everywhere.
#
#   2. TRIAL MIXING. B2-B6 keep NOR habituation and novel-object trials in one
#      folder. Only _NOV files are staged.
#
#   3. DUPLICATE SocP COPIES. B3 and B5 additionally hold 57 flat copies of
#      their HAB/S1/S2 contents. Only the phase subfolders are staged.
#
#   4. EPM `neck`. Added as an alias of spine1, as for Batch 1.
#
# Layout written, matching each script's own path construction:
#   EPM   STAGE/EPM/<batch>/<CODE>.csv
#   NOR   STAGE/NOR/<batch>/<CODE>_NOV.csv   + STAGE/NOR_meta/<batch>/novelLoc.txt
#   SocP  STAGE/SocP/<batch>/<S1|S2>/<CODE>_<phase>.csv + novelLoc<phase>.txt
#   OFT   STAGE/OFT/<batch>/<CODE>.csv
# ============================================================================

suppressMessages({
  library(dplyr)
  library(data.table)
})

RAW   <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior"
PLAN  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Planning"
PROJ  <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
STAGE <- "C:/Users/topohl/Documents/exp9_sleap_staging"
RES   <- file.path(PROJ, "results")
BATCHES <- paste0("B", 1:6)

dir.create(STAGE, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

id_code <- read.delim(file.path(PLAN, "animalIDCode.txt"), stringsAsFactors = FALSE) %>%
  mutate(across(everything(), trimws))
CODES <- id_code$Code

# Known misnamings that no general rule should paper over silently.
ALIASES <- c(M8P1 = "M8P3")

# Resolve the animal a file belongs to. Direct substring match first; then a
# lookalike-character pass (lowercase l for capital I, digit 0 for letter O),
# which is how F2l1/S7l3/R803 are misfiled in the raw BORIS exports; then the
# explicit alias table. Anything still unresolved is reported, never guessed.
resolve_code <- function(stem) {
  hit <- CODES[vapply(CODES, function(cd) grepl(cd, stem, fixed = TRUE), logical(1))]
  if (length(hit) == 1) return(hit)
  if (length(hit) > 1) return(NA_character_)      # ambiguous: report it
  alt <- gsub("l", "I", gsub("0", "O", stem))
  hit <- CODES[vapply(CODES, function(cd) grepl(cd, alt, fixed = TRUE), logical(1))]
  if (length(hit) == 1) return(hit)
  tok <- regmatches(stem, gregexpr("[A-Za-z0-9]{4}", stem))[[1]]
  a <- ALIASES[tok[tok %in% names(ALIASES)]]
  if (length(a) == 1) return(unname(a))
  NA_character_
}

# --- EPM needs a `neck` column added while copying -------------------------
copy_with_neck <- function(src, dst, source_point = "spine1", new_point = "neck") {
  headers <- readLines(src, n = 3L)
  parts <- strsplit(headers, ",", fixed = TRUE)
  bodyparts <- trimws(parts[[2]]); coords <- tolower(trimws(parts[[3]]))
  if (new_point %in% bodyparts) { file.copy(src, dst, overwrite = TRUE); return("copied (had neck)") }
  idx <- which(bodyparts == source_point)
  if (length(idx) != 3L || !identical(coords[idx], c("x", "y", "likelihood")))
    return(sprintf("SKIPPED: no usable %s triplet", source_point))
  d <- fread(src, skip = 3L, header = FALSE, showProgress = FALSE)
  if (ncol(d) != length(bodyparts)) return("SKIPPED: header/data mismatch")
  out <- cbind(d, d[, idx, with = FALSE])
  writeLines(c(
    paste(c(parts[[1]], paste0("column_", new_point, c("_x", "_y", "_lh"))), collapse = ","),
    paste(c(bodyparts, rep(new_point, 3L)), collapse = ","),
    paste(c(parts[[3]], c("x", "y", "likelihood")), collapse = ",")), dst)
  fwrite(out, dst, append = TRUE, col.names = FALSE, showProgress = FALSE)
  "copied (+neck)"
}

log <- list()
note <- function(...) log[[length(log) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

stage_set <- function(assay, batch, src_dir, dst_dir, pattern = "[.]csv$",
                      exclude = NULL, suffix = "", add_neck = FALSE) {
  if (!dir.exists(src_dir)) { note(assay = assay, batch = batch, phase = suffix,
                                   staged = 0L, skipped = 0L, issue = "source missing"); return(invisible()) }
  f <- list.files(src_dir, pattern = pattern, full.names = TRUE)
  if (!is.null(exclude)) f <- f[!grepl(exclude, basename(f), ignore.case = TRUE)]
  dir.create(dst_dir, showWarnings = FALSE, recursive = TRUE)
  staged <- 0L; skipped <- character()
  seen <- character()
  for (p in f) {
    stem <- tools::file_path_sans_ext(basename(p))
    cd <- resolve_code(stem)
    if (is.na(cd)) { skipped <- c(skipped, basename(p)); next }
    if (cd %in% seen) { skipped <- c(skipped, paste0(basename(p), " [dup ", cd, "]")); next }
    seen <- c(seen, cd)
    dst <- file.path(dst_dir, paste0(cd, suffix, ".csv"))
    ok <- if (add_neck) copy_with_neck(p, dst) else
      if (file.copy(p, dst, overwrite = TRUE)) "copied" else "copy failed"
    if (grepl("^copied", ok)) staged <- staged + 1L else skipped <- c(skipped, paste0(basename(p), " [", ok, "]"))
  }
  note(assay = assay, batch = batch, phase = suffix, staged = staged,
       skipped = length(skipped),
       issue = if (length(skipped)) paste(head(skipped, 3), collapse = "; ") else "")
  invisible()
}

cat("staging to:", STAGE, "\n\n")

for (b in BATCHES) {
  cat("---", b, "---\n")

  # EPM: one file per animal, plus the neck alias.
  stage_set("EPM", b, file.path(RAW, b, "EPM/SLEAP/formatted"),
            file.path(STAGE, "EPM", b), add_neck = TRUE)

  # NOR: novel-object trials only, selected by EXCLUDING habituation rather
  # than by requiring "_NOV". Five B5 animals (F8O9, H1E2, M2L3, M7X3, Q3B3)
  # name their novel trial "_NOR" instead of "_NOV"; each has exactly one _HAB
  # and one other file, so the non-HAB file is the novel trial in every batch.
  # Requiring "_NOV" silently dropped those five.
  stage_set("NOR", b, file.path(RAW, b, "NOR/SLEAP/formatted"),
            file.path(STAGE, "NOR", b),
            exclude = "_HAB|combined_output", suffix = "_NOV")
  dir.create(file.path(STAGE, "NOR_meta", b), showWarnings = FALSE, recursive = TRUE)
  file.copy(file.path(RAW, b, "NOR/novelLoc.txt"),
            file.path(STAGE, "NOR_meta", b, "novelLoc.txt"), overwrite = TRUE)

  # SocP: phase subfolders only, never the flat duplicates in B3/B5.
  for (ph in c("S1", "S2")) {
    stage_set("SocP", b, file.path(RAW, b, "SocP/SLEAP/formatted", ph),
              file.path(STAGE, "SocP", b, ph), suffix = paste0("_", ph))
    file.copy(file.path(RAW, b, "SocP", paste0("novelLoc", ph, ".txt")),
              file.path(STAGE, "SocP", b, ph, paste0("novelLoc", ph, ".txt")),
              overwrite = TRUE)
  }

  # OFT: one file per animal.
  stage_set("OFT", b, file.path(RAW, b, "OFT/SLEAP/formatted"),
            file.path(STAGE, "OFT", b))
}

report <- bind_rows(log)
write.csv(report, file.path(RES, "staging_report.csv"), row.names = FALSE)

cat("\n=== STAGED FILE COUNTS ===\n")
print(report %>% mutate(key = paste0(assay, ifelse(nzchar(phase), phase, ""))) %>%
        select(key, batch, staged) %>%
        tidyr::pivot_wider(names_from = batch, values_from = staged, values_fill = 0) %>%
        as.data.frame(), row.names = FALSE)
cat(sprintf("\ntotal staged: %d files\n", sum(report$staged)))

prob <- report %>% filter(skipped > 0 | nzchar(issue))
cat("\n=== SKIPPED / UNRESOLVED ===\n")
if (nrow(prob) == 0) cat("none\n") else print(as.data.frame(prob), row.names = FALSE)

# Shared key file, used by every config. This is the FULL 119-animal map from
# Planning/, not the 20-animal Batch-1 subset in the project's metadata/.
write.table(id_code, file.path(STAGE, "animalIDCode.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nanimalIDCode.txt staged with %d animals\n", nrow(id_code)))
cat("\nwrote:", file.path(RES, "staging_report.csv"), "\n")
