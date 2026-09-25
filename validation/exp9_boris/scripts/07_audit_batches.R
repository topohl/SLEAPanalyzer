# ============================================================================
# 07_audit_batches.R
#
# Pre-flight audit before extending the pipeline from Batch 1 to B1-B6.
# Reads nothing but geometry and file names; runs no analysis and writes
# nothing outside results/.
#
# It exists because the Batch-1 conventions do NOT hold across batches, and
# every way they break is silent:
#
#   1. ANIMAL CODE POSITION. DLCA_*.R takes the code from the first four
#      characters (str_extract(name, "^[A-Za-z0-9]{4}")). That holds for B1,
#      B5 and B6 but not for B2/B3/B4, whose EPM files begin "E9_B2_...",
#      "E9_SIS_B3_..." -- the regex returns NA and the animal silently loses
#      its metadata.
#
#   2. PHASE MIXING. B1's NOR formatted/ folder holds NOV trials only. B2-B6
#      hold HAB and NOV together, so pointing the NOR config at the folder
#      would analyse habituation trials as novel-object trials.
#
#   3. GEOMETRY. Calibration is per-camera. If the rig moved between batches,
#      a single calibration_distance_cm is wrong for some of them.
#
# Reported, not fixed: the fix belongs in the run configuration, once the
# numbers below have been looked at.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
})

PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
RAW  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior"
RES  <- file.path(PROJ, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

id_code <- read.delim("s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Planning/animalIDCode.txt",
                      stringsAsFactors = FALSE) %>% mutate(across(everything(), trimws))
CODES <- id_code$Code
BATCHES <- paste0("B", 1:6)

# The pipeline's own extraction rule, reproduced exactly.
pipeline_code <- function(x) {
  out <- rep(NA_character_, length(x))
  hit <- regexpr("^[A-Za-z0-9]{4}", x)
  out[hit > 0] <- regmatches(x, hit)
  out
}

# What the code actually is: the known code appearing anywhere in the name.
true_code <- function(x) {
  vapply(x, function(nm) {
    found <- CODES[vapply(CODES, function(cd) grepl(cd, nm, fixed = TRUE), logical(1))]
    if (length(found) == 1) found else if (length(found) > 1) paste(found, collapse = "|") else NA_character_
  }, character(1), USE.NAMES = FALSE)
}

# --- 1. Inventory and code extraction --------------------------------------
assay_dirs <- function(batch, assay) {
  base <- file.path(RAW, batch, assay, "SLEAP", "formatted")
  if (!dir.exists(base)) return(character(0))
  subs <- list.dirs(base, recursive = FALSE)
  subs <- subs[!grepl("/old$", subs)]
  if (length(subs) > 0) c(base, subs) else base
}

inv <- list()
for (b in BATCHES) for (a in c("EPM", "NOR", "SocP", "OFT")) {
  for (d in assay_dirs(b, a)) {
    f <- list.files(d, pattern = "[.]csv$")
    if (length(f) == 0) next
    phase <- basename(d)
    inv[[length(inv) + 1]] <- data.frame(
      batch = b, assay = a,
      folder = if (phase == "formatted") "(root)" else phase,
      file = f, stringsAsFactors = FALSE)
  }
}
inv <- bind_rows(inv) %>%
  mutate(stem = tools::file_path_sans_ext(file),
         code_pipeline = pipeline_code(stem),
         code_true = true_code(stem),
         code_ok = !is.na(code_pipeline) & !is.na(code_true) &
                   code_pipeline == code_true,
         trial = case_when(grepl("_HAB", file, ignore.case = TRUE) ~ "HAB",
                           grepl("_NOV", file, ignore.case = TRUE) ~ "NOV",
                           TRUE ~ "-"))

write.csv(inv, file.path(RES, "batch_audit_files.csv"), row.names = FALSE)

cat("=== 1. FILE INVENTORY ===\n")
print(inv %>% count(assay, batch, folder, name = "files") %>%
        pivot_wider(names_from = batch, values_from = files, values_fill = 0) %>%
        as.data.frame(), row.names = FALSE)

cat("\n=== 2. ANIMAL CODE EXTRACTION (pipeline rule vs reality) ===\n")
codesum <- inv %>%
  group_by(assay, batch) %>%
  summarise(files = n(),
            pipeline_ok = sum(code_ok),
            pipeline_NA = sum(is.na(code_pipeline)),
            wrong_code  = sum(!is.na(code_pipeline) & !code_ok),
            code_unknown = sum(is.na(code_true)), .groups = "drop")
print(as.data.frame(codesum), row.names = FALSE)
cat("\npipeline_ok = the first four characters really are the animal code.\n")
cat("pipeline_NA / wrong_code = the run would lose or mislabel those animals.\n")

bad <- codesum %>% filter(pipeline_ok < files)
if (nrow(bad)) {
  cat("\n--- example broken names ---\n")
  ex <- inv %>% filter(!code_ok) %>% group_by(assay, batch) %>% slice(1) %>% ungroup()
  print(ex %>% select(assay, batch, file, code_pipeline, code_true) %>%
          as.data.frame(), row.names = FALSE)
}

unk <- inv %>% filter(is.na(code_true))
if (nrow(unk)) {
  cat("\n--- files whose animal could not be identified at all ---\n")
  print(unk %>% select(assay, batch, folder, file) %>% as.data.frame(), row.names = FALSE)
}

cat("\n=== 3. TRIAL MIXING (NOR) ===\n")
print(inv %>% filter(assay == "NOR") %>% count(batch, folder, trial) %>%
        pivot_wider(names_from = trial, values_from = n, values_fill = 0) %>%
        as.data.frame(), row.names = FALSE)
cat("A batch with HAB and NOV in the same folder needs the config restricted\n",
    "to NOV, or habituation trials are scored as novel-object trials.\n", sep = "")

cat("\n=== 4. DUPLICATE ANIMALS WITHIN A FOLDER ===\n")
dup <- inv %>% filter(!is.na(code_true)) %>%
  count(assay, batch, folder, trial, code_true) %>% filter(n > 1)
if (nrow(dup) == 0) cat("none\n") else print(as.data.frame(dup), row.names = FALSE)

cat("\n=== 5. ANIMALS PER BATCH vs the metadata tables ===\n")
# Stage 01 foundation copy, identical to the numbered original being archived.
asg <- read.csv(file.path("s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress",
                          "Analysis/Behavior/RFID/analysis_ready/foundations/behavior_metrics/qc",
                          "animal_group_sex_assignment_qc.csv"), stringsAsFactors = FALSE)
cov <- inv %>% filter(!is.na(code_true), assay %in% c("EPM", "NOR", "SocP")) %>%
  distinct(batch, assay, code_true) %>%
  left_join(id_code, by = c("code_true" = "Code")) %>%
  mutate(in_rfid = ID %in% asg$AnimalID_norm |
                   sub("^O[QR]0*", "", ID) %in% asg$AnimalID_norm)
print(cov %>% group_by(assay, batch) %>%
        summarise(animals = n(), in_animalIDCode = sum(!is.na(ID)),
                  in_rfid_table = sum(in_rfid), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\nwrote:", file.path(RES, "batch_audit_files.csv"), "\n")

# --- 6. Geometry / calibration stability across batches --------------------
# Calibration is per-camera. If the rig moved between batches, one
# calibration_distance_cm cannot serve all of them. Only the static geometry
# landmarks are needed, so a few hundred frames per file is plenty.
suppressMessages(library(data.table))

geom_medians <- function(path, want, nrows = 300L) {
  hdr <- readLines(path, n = 3L)
  parts <- strsplit(hdr, ",", fixed = TRUE)
  bodyparts <- trimws(parts[[2]]); coords <- tolower(trimws(parts[[3]]))
  d <- tryCatch(fread(path, skip = 3L, header = FALSE, nrows = nrows,
                      showProgress = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  out <- lapply(want, function(p) {
    ix <- which(bodyparts == p & coords %in% c("x", "y"))
    if (length(ix) < 2) return(c(x = NA_real_, y = NA_real_))
    c(x = median(as.numeric(d[[ix[1]]]), na.rm = TRUE),
      y = median(as.numeric(d[[ix[2]]]), na.rm = TRUE))
  })
  names(out) <- want
  out
}
dist2 <- function(a, b) sqrt((a[["x"]] - b[["x"]])^2 + (a[["y"]] - b[["y"]])^2)

cat("\n=== 6. GEOMETRY BY BATCH ===\n")

cat("\n--- EPM: tl-bl is the 60 cm tip-to-tip span; arm width should read ~5 cm ---\n")
epm_geo <- bind_rows(lapply(BATCHES, function(b) {
  f <- list.files(file.path(RAW, b, "EPM/SLEAP/formatted"),
                  pattern = "[.]csv$", full.names = TRUE)
  f <- head(f, 3)
  bind_rows(lapply(f, function(p) {
    g <- geom_medians(p, c("tl", "tr", "bl", "br"))
    if (is.null(g) || any(is.na(g$tl))) return(NULL)
    s <- 60 / dist2(g$tl, g$bl)
    data.frame(batch = b, file = basename(p),
               tl_bl_px = dist2(g$tl, g$bl), tl_br_px = dist2(g$tl, g$br),
               cm_per_px = s, arm_width_cm = dist2(g$tl, g$tr) * s,
               stringsAsFactors = FALSE)
  }))
}))
print(epm_geo %>% group_by(batch) %>%
        summarise(n = n(), tl_bl_px = round(median(tl_bl_px), 1),
                  cm_per_px = round(median(cm_per_px), 5),
                  arm_width_cm = round(median(arm_width_cm), 2), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\n--- NOR / SocP: arena corner polygon, as a rectangle check ---\n")
rect_geo <- bind_rows(lapply(c("NOR", "SocP"), function(a) {
  bind_rows(lapply(BATCHES, function(b) {
    base <- file.path(RAW, b, a, "SLEAP/formatted")
    d <- if (a == "SocP" && dir.exists(file.path(base, "S1"))) file.path(base, "S1") else base
    f <- head(list.files(d, pattern = "[.]csv$", full.names = TRUE), 3)
    bind_rows(lapply(f, function(p) {
      g <- geom_medians(p, c("tl", "tr", "br", "bl"))
      if (is.null(g) || any(is.na(g$tl))) return(NULL)
      w <- mean(c(dist2(g$tl, g$tr), dist2(g$bl, g$br)))
      h <- mean(c(dist2(g$tl, g$bl), dist2(g$tr, g$br)))
      data.frame(assay = a, batch = b, width_px = w, height_px = h,
                 aspect = w / h, stringsAsFactors = FALSE)
    }))
  }))
}))
print(rect_geo %>% group_by(assay, batch) %>%
        summarise(n = n(), width_px = round(median(width_px), 1),
                  height_px = round(median(height_px), 1),
                  aspect = round(median(aspect), 3), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("\nNOR arena is square, so aspect should be ~1.00; SocP is 44 x 24 cm, ~1.83.\n")
cat("A batch whose px sizes differ markedly from the others had the camera moved\n",
    "and needs its own calibration values.\n", sep = "")

write.csv(epm_geo, file.path(RES, "batch_audit_geometry_epm.csv"), row.names = FALSE)
write.csv(rect_geo, file.path(RES, "batch_audit_geometry_rect.csv"), row.names = FALSE)
