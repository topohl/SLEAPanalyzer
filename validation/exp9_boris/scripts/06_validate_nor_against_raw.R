# ============================================================================
# 06_validate_nor_against_raw.R
#
# NOR.xlsx is a transcription. The primary records are the 40 per-animal BORIS
# exports in Raw Data/Behavior/B1/NOR/BORIS/ (*_hab.tsv, *_nov.tsv), which
# carry "Interaction left" / "Interaction Right" straight out of BORIS.
#
# This script checks the sheet against those exports. Two things came out of
# it, and both change how the NOR results should be read:
#
#   1. The sheet's intLeft_BORIS / intRight_BORIS are SWAPPED relative to the
#      raw export's own left/right labels, for every animal that reconciles.
#      That is the source of the "mirror" seen against SLEAP: measured against
#      the RAW export, SLEAP's contactLeft agrees with "Interaction left"
#      directly (r = +0.95), so neither SLEAP nor the scorer was mirrored --
#      the sheet's two column headers are.
#
#   2. Three sheet rows (A1L9, F9L3, R8K7) match no value in any raw export,
#      so they came from a different scoring pass. They alone drag D2
#      agreement from r = 0.960 down to r = 0.740.
#
# All 20 animals have a *_nov.tsv export. Three of them are misnamed with
# lookalike characters, so files are matched on their Subject column rather
# than on the file name -- see the note above the reader below.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(readxl)
})

PROJ  <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
BORIS <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B1/NOR/BORIS"
RES   <- file.path(PROJ, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

# --- Read the raw per-animal exports ---------------------------------------
# Key on the file's own Subject column, NOT the file name. Three exports are
# misnamed with lookalike characters -- F2l1_nov.tsv and S7l3_nov.tsv use a
# lowercase L for the capital I, and R803_nov.tsv a digit 0 for the letter O
# (the same slip appears in the NOR .slp video names). The Subject column
# inside each file is correct, so matching on it recovers all 20 animals;
# matching on the file name silently loses those three.
files <- list.files(BORIS, pattern = "_nov[.]tsv$", full.names = TRUE)
stopifnot(length(files) > 0)
raw <- do.call(rbind, lapply(files, function(p) {
  d <- read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
  pick <- function(cat) {
    v <- d[["Total duration (s)"]][trimws(d$Category) == cat]
    if (length(v) == 0) NA_real_ else v[1]
  }
  n_of <- function(cat) {
    v <- d[["Total number of occurences"]][trimws(d$Category) == cat]
    if (length(v) == 0) NA_integer_ else as.integer(v[1])
  }
  subject <- unique(trimws(d$Subject))
  subject <- subject[nzchar(subject)]
  if (length(subject) != 1) {
    stop(basename(p), ": expected one Subject, found ",
         paste(subject, collapse = ", "))
  }
  data.frame(
    Code = subject,
    file_name = basename(p),
    name_matches_subject = identical(subject, sub("_nov$", "", tools::file_path_sans_ext(basename(p)))),
    raw_left = pick("Interaction left"), raw_right = pick("Interaction Right"),
    raw_left_n = n_of("Interaction left"), raw_right_n = n_of("Interaction Right"),
    stringsAsFactors = FALSE
  )
}))
stopifnot(!any(duplicated(raw$Code)))
if (any(!raw$name_matches_subject)) {
  cat("misnamed raw exports (Subject column used instead):\n")
  print(raw %>% filter(!name_matches_subject) %>% select(file_name, Code),
        row.names = FALSE)
  cat("\n")
}

sheet <- suppressMessages(read_excel(file.path(PROJ, "NOR.xlsx"))) %>%
  transmute(Code = trimws(Subject), novelloc = trimws(novelloc),
            intLeft_BORIS, intRight_BORIS, nov_BORIS, fam_BORIS, D2_BORIS)

# --- Reconciliation --------------------------------------------------------
near <- function(a, b) !is.na(a) & !is.na(b) & abs(a - b) < 1e-3

recon <- sheet %>%
  left_join(raw, by = "Code") %>%
  mutate(
    has_raw       = !is.na(raw_left),
    sides_direct  = near(intLeft_BORIS, raw_left)  & near(intRight_BORIS, raw_right),
    sides_swapped = near(intLeft_BORIS, raw_right) & near(intRight_BORIS, raw_left),
    novfam_from_raw =
      (near(nov_BORIS, raw_left)  & near(fam_BORIS, raw_right)) |
      (near(nov_BORIS, raw_right) & near(fam_BORIS, raw_left)),
    # Which raw (image) side did the sheet treat as novel?
    nov_raw_side = case_when(near(nov_BORIS, raw_left)  ~ "L",
                             near(nov_BORIS, raw_right) ~ "R",
                             TRUE ~ NA_character_),
    status = case_when(
      !has_raw                      ~ "no raw export",
      sides_swapped & novfam_from_raw ~ "reconciles (sheet sides swapped)",
      sides_direct  & novfam_from_raw ~ "reconciles (sheet sides direct)",
      TRUE                          ~ "DOES NOT RECONCILE"
    )
  )

write.csv(recon, file.path(RES, "nor_sheet_vs_raw_reconciliation.csv"), row.names = FALSE)

cat("=== NOR.xlsx vs raw BORIS exports ===\n")
print(recon %>% count(status) %>% as.data.frame(), row.names = FALSE)

cat("\n--- sheet column orientation, among animals with a raw export ---\n")
cat(sprintf("  intLeft == raw 'Interaction left'  : %d\n", sum(recon$sides_direct, na.rm = TRUE)))
cat(sprintf("  intLeft == raw 'Interaction Right' : %d  <- the swap\n",
            sum(recon$sides_swapped, na.rm = TRUE)))

cat("\n--- novel side, relative to the RAW (image) sides ---\n")
print(recon %>% filter(!is.na(nov_raw_side)) %>%
        count(novelloc, nov_raw_side) %>% as.data.frame(), row.names = FALSE)
cat("novelloc 'R' -> novel is raw LEFT: the same convention the pipeline uses\n",
    "(novel_is_left <- location == 'R'), which is why contactNov matches nov_BORIS.\n", sep = "")

cat("\n--- rows needing attention ---\n")
print(recon %>% filter(status %in% c("DOES NOT RECONCILE", "no raw export")) %>%
        select(Code, status, novelloc, raw_left, raw_right,
               intLeft_BORIS, intRight_BORIS, nov_BORIS, fam_BORIS) %>%
        as.data.frame(), row.names = FALSE)

# --- Agreement against raw, and the cost of the three bad rows -------------
sleap <- read.csv(file.path(PROJ, "sleap_output/NOR/combined_output.csv"),
                  stringsAsFactors = FALSE) %>%
  mutate(Code = substr(basename(file), 1, 4)) %>%
  select(Code, cL = contactLeft, cR = contactRight,
         cNov = contactNov, cFam = contactFam, novelLoc)

ccc <- function(x, y) {
  n <- length(x); vx <- var(x) * (n - 1) / n; vy <- var(y) * (n - 1) / n
  2 * cov(x, y) * (n - 1) / n / (vx + vy + (mean(x) - mean(y))^2)
}

cmp <- recon %>%
  filter(has_raw) %>%
  inner_join(sleap, by = "Code") %>%
  mutate(
    raw_nov  = ifelse(novelLoc == "R", raw_left, raw_right),
    raw_fam  = ifelse(novelLoc == "R", raw_right, raw_left),
    raw_D2   = (raw_nov - raw_fam) / (raw_nov + raw_fam),
    sleap_D2 = (cNov - cFam) / (cNov + cFam)
  )

row_of <- function(d, label, sleap_v, ref_v, ref) {
  data.frame(comparison = label, reference = ref, n = nrow(d),
             r = cor(sleap_v, ref_v), ccc = ccc(sleap_v, ref_v),
             bias = mean(sleap_v - ref_v), stringsAsFactors = FALSE)
}
keep <- cmp %>% filter(status != "DOES NOT RECONCILE")

summary_tbl <- bind_rows(
  row_of(cmp,  "contactLeft ~ Interaction left",  cmp$cL,  cmp$raw_left,  "raw export"),
  row_of(cmp,  "contactRight ~ Interaction Right",cmp$cR,  cmp$raw_right, "raw export"),
  row_of(cmp,  "contactNov ~ novel",              cmp$cNov, cmp$raw_nov,  "raw export"),
  row_of(cmp,  "D2",                              cmp$sleap_D2, cmp$raw_D2, "raw export"),
  row_of(cmp,  "contactNov ~ novel",              cmp$cNov, cmp$nov_BORIS, "NOR.xlsx"),
  row_of(cmp,  "D2",                              cmp$sleap_D2, cmp$D2_BORIS, "NOR.xlsx"),
  row_of(keep, "contactNov ~ novel",              keep$cNov, keep$raw_nov, "raw export, reconciling only"),
  row_of(keep, "D2",                              keep$sleap_D2, keep$raw_D2, "raw export, reconciling only")
)
write.csv(summary_tbl, file.path(RES, "nor_agreement_vs_raw.csv"), row.names = FALSE)

cat("\n=== SLEAP agreement, by reference ===\n")
print(summary_tbl %>% mutate(across(c(r, ccc), ~round(., 3)), bias = round(bias, 2)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nSLEAP's contactLeft agrees with the raw 'Interaction left' directly, so the\n",
    "mirror lives in the NOR.xlsx headers, not in SLEAP or in the scoring.\n", sep = "")
cat("Dropping the three non-reconciling rows raises D2 agreement from r = ",
    sprintf("%.3f to %.3f.\n", summary_tbl$r[summary_tbl$comparison == "D2" &
                                             summary_tbl$reference == "raw export"],
            summary_tbl$r[summary_tbl$reference == "raw export, reconciling only" &
                          summary_tbl$comparison == "D2"]), sep = "")
