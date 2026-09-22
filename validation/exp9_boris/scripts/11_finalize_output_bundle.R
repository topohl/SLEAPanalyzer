# ============================================================================
# 11_finalize_output_bundle.R
#
# Finalise an immutable, canonical-only Exp9 SLEAPanalyzer release. The
# assembly and analysis scripts must already have written into the same run
# root. This script creates clearly named all-animal and sex-specific handoff
# tables plus compact coverage/missingness QA and a human-readable README.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(writexl)
})

RUN_ROOT <- Sys.getenv("EXP9_SLEAP_RUN_ROOT")
if (!nzchar(RUN_ROOT)) {
  stop("EXP9_SLEAP_RUN_ROOT must name the release directory.")
}
DATA_DIR <- Sys.getenv("EXP9_SLEAP_DATA_DIR", unset = file.path(RUN_ROOT, "data"))
META_DIR <- Sys.getenv("EXP9_SLEAP_METADATA_DIR", unset = file.path(RUN_ROOT, "metadata"))
QC_DIR <- file.path(RUN_ROOT, "qc")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(META_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

assembled_tsv <- file.path(DATA_DIR, "sleap_all_batches_wide.tsv")
assembled_xlsx <- file.path(DATA_DIR, "sleap_all_batches_wide.xlsx")
metadata_tsv <- file.path(META_DIR, "animal_metadata_all.tsv")
for (path in c(assembled_tsv, assembled_xlsx, metadata_tsv)) {
  if (!file.exists(path)) stop("Required assembly output is missing: ", path)
}

dat <- read.delim(assembled_tsv, stringsAsFactors = FALSE, na.strings = "")
required <- c("Code", "ID", "Batch", "Sex", "Condition", "Phenotype")
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) {
  stop("Assembled table is missing required columns: ",
       paste(missing_required, collapse = ", "))
}
forbidden <- c("Phenotype_batchCorrected", "Phenotype_conflict")
present_forbidden <- intersect(forbidden, names(dat))
if (length(present_forbidden)) {
  stop("Canonical bundle must not contain deprecated phenotype columns: ",
       paste(present_forbidden, collapse = ", "))
}
if (nrow(dat) != 117L || anyDuplicated(dat$Code)) {
  stop("Expected 117 unique animal codes; found ", nrow(dat),
       " rows and ", anyDuplicated(dat$Code), " duplicated-code indicator.")
}

sex_batches <- list(
  Male = c("B1", "B2", "B5"),
  Female = c("B3", "B4", "B6")
)
for (sex in names(sex_batches)) {
  observed <- sort(unique(dat$Batch[dat$Sex == sex]))
  expected <- sort(sex_batches[[sex]])
  if (!identical(observed, expected)) {
    stop(sex, " batches are ", paste(observed, collapse = ", "),
         "; expected ", paste(expected, collapse = ", "))
  }
}
if (any(is.na(dat$Sex)) || any(is.na(dat$Condition)) || any(is.na(dat$Phenotype))) {
  stop("Sex, Condition and Phenotype must be complete for all bundled animals.")
}

male <- dat %>% filter(Sex == "Male")
female <- dat %>% filter(Sex == "Female")
stopifnot(nrow(male) == 59L, nrow(female) == 58L)
expected_condition <- c(CON = 24L, SIS = 93L)
expected_phenotype <- c(CON = 24L, RES = 58L, SUS = 35L)
observed_condition <- table(dat$Condition)
observed_phenotype <- table(dat$Phenotype)
if (!setequal(names(observed_condition), names(expected_condition)) ||
    !identical(as.integer(observed_condition[names(expected_condition)]),
               unname(expected_condition))) {
  stop("Unexpected condition counts: ",
       paste(names(observed_condition), observed_condition, collapse = ", "))
}
if (!setequal(names(observed_phenotype), names(expected_phenotype)) ||
    !identical(as.integer(observed_phenotype[names(expected_phenotype)]),
               unname(expected_phenotype))) {
  stop("Unexpected canonical phenotype counts: ",
       paste(names(observed_phenotype), observed_phenotype, collapse = ", "))
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

all_path <- file.path(DATA_DIR, "exp9_sleap_all_animals.tsv")
male_path <- file.path(DATA_DIR, "exp9_sleap_males_B1-B2-B5.tsv")
female_path <- file.path(DATA_DIR, "exp9_sleap_females_B3-B4-B6.tsv")
meta_path <- file.path(META_DIR, "exp9_sleap_animal_metadata.tsv")
xlsx_path <- file.path(DATA_DIR, "exp9_sleap_analysis_tables.xlsx")

write_tsv(dat, all_path)
write_tsv(male, male_path)
write_tsv(female, female_path)
if (!file.copy(metadata_tsv, meta_path, overwrite = FALSE)) {
  stop("Could not create canonical metadata handoff: ", meta_path)
}
write_xlsx(
  list(All_animals = dat, Male_B1_B2_B5 = male, Female_B3_B4_B6 = female),
  xlsx_path
)

assay_columns <- c(
  EPM = "sleap_EPM_open_s",
  NOR = "sleap_NOR_nov_s",
  SocP = "sleap_SocP_S1_novel",
  OFT = "sleap_OFT_center_s"
)
coverage_long <- bind_rows(lapply(names(assay_columns), function(assay) {
  column <- assay_columns[[assay]]
  if (!column %in% names(dat)) stop("Missing assay coverage column: ", column)
  data.frame(
    Code = dat$Code,
    ID = dat$ID,
    Batch = dat$Batch,
    Sex = dat$Sex,
    assay = assay,
    available = !is.na(dat[[column]]),
    stringsAsFactors = FALSE
  )
}))
coverage <- coverage_long %>%
  group_by(assay, Batch, Sex) %>%
  summarise(animals = n(), available = sum(available),
            missing = sum(!available), .groups = "drop")
missing <- coverage_long %>% filter(!available) %>% select(-available)
cohort <- dat %>%
  count(Batch, Sex, Condition, Phenotype, name = "animals") %>%
  arrange(Sex, Batch, Condition, Phenotype)

write.csv(coverage, file.path(QC_DIR, "assay_coverage_by_batch.csv"), row.names = FALSE)
write.csv(missing, file.path(QC_DIR, "missing_assay_records.csv"), row.names = FALSE)
write.csv(cohort, file.path(QC_DIR, "cohort_counts.csv"), row.names = FALSE)

readme <- c(
  "# Exp9 SLEAPanalyzer v2 all-batch release",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Cohort",
  "",
  "- 117 unique animals across B1-B6.",
  "- Male: 59 animals from B1, B2 and B5.",
  "- Female: 58 animals from B3, B4 and B6.",
  "- Condition: 24 CON and 93 SIS.",
  sprintf("- Canonical phenotype: %d CON, %d RES and %d SUS.",
          sum(dat$Phenotype == "CON"), sum(dat$Phenotype == "RES"),
          sum(dat$Phenotype == "SUS")),
  "- Deprecated batch-corrected phenotype lists are not read or bundled.",
  "",
  "Sex is perfectly nested in batch. Sex-stratified models combine three",
  "batches within each sex and retain Batch as a blocking factor. A sex main",
  "effect adjusted for six-level Batch is not estimable.",
  "",
  "## Directory map",
  "",
  "- data/: canonical all-animal table, explicit male/female tables, workbook.",
  "- metadata/: animal identity and experimental-group metadata.",
  "- assay_summaries/: compact per-batch assay tables, QC and run manifests.",
  "- statistics/: pooled and sex-stratified effect tables/workbook.",
  "- figures/: figures generated from this release's assembled table.",
  "- qc/: cohort counts, assay coverage and missing-record register.",
  "- provenance/: exact configs/scripts, source-copy map, Git/R records and hashes.",
  "",
  "Bulk formatted coordinates, per-animal plots and EPM TIFF overview images",
  "are intentionally excluded. They remain in the recorded source locations.",
  "",
  "## Primary handoff files",
  "",
  "- data/exp9_sleap_all_animals.tsv",
  "- data/exp9_sleap_males_B1-B2-B5.tsv",
  "- data/exp9_sleap_females_B3-B4-B6.tsv",
  "- data/exp9_sleap_analysis_tables.xlsx",
  "- statistics/all_batches_effects.csv",
  "- statistics/all_batches_primary.csv",
  "- statistics/all_batches_sex_interaction.csv",
  "- statistics/all_batches_results.xlsx",
  "- statistics/exp9_sleap_effects_all_animals_batch_adjusted.csv",
  "- statistics/exp9_sleap_effects_males_B1-B2-B5_batch_adjusted.csv",
  "- statistics/exp9_sleap_effects_females_B3-B4-B6_batch_adjusted.csv"
)
writeLines(readme, file.path(RUN_ROOT, "README.md"))

cat(sprintf(
  "bundle finalised: %d animals (%d male, %d female); missing assay records: %d\n",
  nrow(dat), nrow(male), nrow(female), nrow(missing)
))
