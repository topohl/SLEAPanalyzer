# ============================================================================
# 03_prepare_sleap_inputs.R
#
# Prepares the staged SLEAP tracking exports in sleap_input/ so SLEAPanalyzer
# v2 can read them. Two fixes, both applied to the local copies only -- the
# files under Raw Data/ are never touched.
#
# 1. Add a `neck` landmark to the EPM files.
#
#    DLCA_EPM v1.0.0.R lists `neck` in required_points and stops if it is
#    absent. The EPM SLEAP skeleton has 15 nodes and no node of that name, so
#    every file would fail. The skeleton's edge list is
#
#        head_center -- spine_1,  spine_1 -- center,
#        spine_1 -- left_ear / right_ear / left_arm / right_arm
#
#    i.e. spine_1 sits between the head and the body centre and carries the
#    ear and forelimb attachments: it *is* the neck/shoulder landmark. The
#    nose-dip test needs the headcentre -> neck -> bodycentre chain, which
#    maps onto head_center -> spine_1 -> center. `neck` is therefore added as
#    an alias of spine1 rather than as a derived midpoint.
#
# 2. Fix the misnamed EPM file M8P1 -> M8P3.
#
#    The pipeline takes the animal code from the first four characters of the
#    file name, so `M8P1...csv` would produce code "M8P1", which matches no
#    animal in the metadata and would silently lose that animal's grouping.
#    M8P3 is the Batch-1 animal (OQ763); no M8P1 exists in animalIDCode.txt.
# ============================================================================

suppressMessages(library(data.table))

PROJ  <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
INPUT <- file.path(PROJ, "sleap_input")

# --- 1. Rename the misnamed EPM file ---------------------------------------
bad <- list.files(file.path(INPUT, "EPM"), pattern = "^M8P1", full.names = TRUE)
for (f in bad) {
  target <- file.path(dirname(f), sub("^M8P1", "M8P3", basename(f)))
  if (file.exists(target)) {
    message("target already exists, leaving in place: ", basename(target))
  } else {
    file.rename(f, target)
    message("renamed: ", basename(f), " -> ", basename(target))
  }
}

# --- 2. Add the `neck` alias to every EPM file -----------------------------
add_neck_alias <- function(path, source_point = "spine1", new_point = "neck") {
  headers <- readLines(path, n = 3L)
  parts   <- strsplit(headers, ",", fixed = TRUE)

  bodyparts <- trimws(parts[[2]])
  coords    <- tolower(trimws(parts[[3]]))

  if (new_point %in% bodyparts) return("already_present")

  # Locate the source triplet. Point names sit at columns 2, 5, 8, ... and the
  # reader requires the x/y/likelihood order, so the triplet is contiguous.
  idx <- which(bodyparts == source_point)
  if (length(idx) != 3L) {
    stop(basename(path), ": expected exactly 3 '", source_point,
         "' columns, found ", length(idx))
  }
  if (!identical(coords[idx], c("x", "y", "likelihood"))) {
    stop(basename(path), ": '", source_point,
         "' columns are not in x/y/likelihood order")
  }

  dat <- fread(path, skip = 3L, header = FALSE, showProgress = FALSE)
  if (ncol(dat) != length(bodyparts)) {
    stop(basename(path), ": header/data column count mismatch")
  }

  # Append the copied triplet, and extend all three header rows to match.
  dat_out <- cbind(dat, dat[, idx, with = FALSE])
  new_scorer    <- c(parts[[1]], paste0("column_", new_point, c("_x", "_y", "_lh")))
  new_bodyparts <- c(bodyparts, rep(new_point, 3L))
  new_coords    <- c(parts[[3]], c("x", "y", "likelihood"))

  tmp <- paste0(path, ".tmp")
  writeLines(c(paste(new_scorer,    collapse = ","),
               paste(new_bodyparts, collapse = ","),
               paste(new_coords,    collapse = ",")), tmp)
  fwrite(dat_out, tmp, append = TRUE, col.names = FALSE, showProgress = FALSE)
  file.rename(tmp, path)
  "added"
}

epm_files <- list.files(file.path(INPUT, "EPM"), pattern = "\\.csv$", full.names = TRUE)
cat(sprintf("EPM files: %d\n", length(epm_files)))
res <- vapply(epm_files, add_neck_alias, character(1))
print(table(res))

# --- 3. Verify every staged file is readable by the pipeline ---------------
source("C:/Users/topohl/Documents/GitHub/SLEAPanalyzer/02_SLEAPanalzyer/DLCAnalyzer_Functions_final.R")

check <- function(path, need) {
  t <- tryCatch(ReadDLCDataFromCSV(path, fps = 30),
                error = function(e) conditionMessage(e))
  if (is.character(t)) return(data.frame(file = basename(path), ok = FALSE,
                                         detail = t, stringsAsFactors = FALSE))
  missing <- setdiff(need, names(t$data))
  data.frame(
    file   = basename(path),
    ok     = length(missing) == 0,
    detail = if (length(missing)) paste("missing:", paste(missing, collapse = ", "))
             else sprintf("%d points, %d frames", length(t$data), length(t$frames)),
    stringsAsFactors = FALSE
  )
}

epm_need  <- c("tl", "tr", "br", "bl", "ctr", "rt", "rb", "cbr", "bl", "cbl",
               "lb", "lt", "ctl", "headcentre", "bodycentre", "neck")
nor_need  <- c("tl", "tr", "br", "bl", "objL", "objR", "nose", "bodycentre")
socp_need <- c("tl", "tr", "br", "bl", "socl", "socr", "nose", "bodycentre")

cat("\n=== readability / landmark check ===\n")
for (set in list(
  list(dir = "EPM",     need = epm_need),
  list(dir = "NOR",     need = nor_need),
  list(dir = "SocP/S1", need = socp_need),
  list(dir = "SocP/S2", need = socp_need)
)) {
  files <- list.files(file.path(INPUT, set$dir), pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) { cat(sprintf("%-8s -- no files staged\n", set$dir)); next }
  out <- do.call(rbind, lapply(files, check, need = set$need))
  cat(sprintf("%-8s %d files | all readable: %s\n",
              set$dir, nrow(out), all(out$ok)))
  if (!all(out$ok)) print(out[!out$ok, ], row.names = FALSE)

  codes <- substr(basename(files), 1, 4)
  meta  <- read.delim(file.path(PROJ, "metadata/animal_metadata.tsv"),
                      stringsAsFactors = FALSE, na.strings = "")
  unmatched <- setdiff(codes, meta$Code)
  cat(sprintf("         codes matching metadata: %d/%d%s\n",
              sum(codes %in% meta$Code), length(codes),
              if (length(unmatched)) paste0(" | UNMATCHED: ", paste(unmatched, collapse = ", ")) else ""))
}
