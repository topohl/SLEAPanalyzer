# ============================================================================
# 01_build_metadata.R
#
# Builds a canonical animal-level metadata table for the Exp9 Batch-1 animals
# that have manual (BORIS) tracking data, and writes metadata-enriched copies
# of the NOR, SocP and EPM scoring files.
#
# Design rules (from Exp9_Social-Stress/README.md):
#   - Original experimental records are never modified; enriched copies only.
#   - Experimental condition (CON/SIS) is kept distinct from derived
#     phenotype (RES/SUS).
#   - Every derived classification carries the source it came from.
#   - Conflicting metadata are recorded, not silently corrected.
#
# Join key is `Code` (the four-character animal code, e.g. A1L9), which is
# also the key SLEAPanalyzer metadata_lookup() uses.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(readxl)
  library(writexl)
  library(tidyr)
})

# --- Paths ------------------------------------------------------------------
EXP9   <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress"
PROJ   <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
META   <- file.path(PROJ, "metadata")
ENRICH <- file.path(PROJ, "enriched")
dir.create(META,   showWarnings = FALSE, recursive = TRUE)
dir.create(ENRICH, showWarnings = FALSE, recursive = TRUE)

# The phenotype table exists only in the retained numbered original
# 03_derived_metrics, which may be archived under history/original_layout/.
# Its archive receipt, not directory existence, says where it is; any other
# state stops the script.
retained_derived_metrics <- function(exp9) {
  ready <- file.path(exp9, "Analysis/Behavior/RFID/analysis_ready")
  original <- file.path(ready, "03_derived_metrics")
  archived <- file.path(ready, "history/original_layout/03_derived_metrics")
  receipt <- file.path(ready, "_migration_control/numbered_root_archive/03_derived_metrics.json")
  if (!file.exists(receipt)) {
    if (dir.exists(archived)) stop("03_derived_metrics archive exists without a receipt: ", archived)
    return(original)
  }
  state <- jsonlite::fromJSON(receipt, simplifyVector = FALSE)[["state"]]
  if (identical(state, "prepared") && dir.exists(original) && !dir.exists(archived)) return(original)
  if (identical(state, "activated") && !dir.exists(original) && dir.exists(archived)) return(archived)
  stop("03_derived_metrics archive is not readable (state ", format(state), "): ", receipt)
}

src <- list(
  id_code   = file.path(EXP9, "Planning/animalIDCode.txt"),
  # Stage 01 foundation copy, identical to the numbered original.
  assign    = file.path(EXP9, "Analysis/Behavior/RFID/analysis_ready/foundations/behavior_metrics/qc/animal_group_sex_assignment_qc.csv"),
  phenotype = file.path(retained_derived_metrics(EXP9), "qc/cross_scale_identity_expected_phenotype_from_preprocessed.csv"),
  sus       = file.path(EXP9, "Analysis/sus_animals.csv"),
  sus_bc    = file.path(EXP9, "Analysis/sus_animals_batchCorrected.csv"),
  con       = file.path(EXP9, "Analysis/con_animals.csv"),
  nor       = file.path(PROJ, "NOR.xlsx"),
  socp      = file.path(PROJ, "SocP.xlsx"),
  epm       = file.path(PROJ, "SISB1_EPMdata.tsv"),
  epm_grp   = file.path(PROJ, "SISB1_EPMdataGrouped.tsv")
)
stopifnot(all(file.exists(unlist(src))))

# --- Identifier normalisation ----------------------------------------------
# Exp9 identifiers appear in several shapes across files because Excel drops
# leading zeros and some exports drop the OQ/OR prefix:
#   animalIDCode.txt : "OQ752", "0003"
#   NOR.xlsx         : "752",   "0003"   (text, prefix dropped)
#   SocP.xlsx        : 752,     3        (numeric, prefix AND padding dropped)
# normalise_id() maps any of these onto the canonical animalIDCode ID by
# trying each plausible reconstruction and keeping the one that exists.
normalise_id <- function(x, canonical) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    v <- trimws(as.character(v))
    cand <- unique(c(
      v,
      sub("\\.0$", "", v),
      if (grepl("^[0-9]+$", v)) c(
        sprintf("%04d", as.integer(v)),   # 3   -> 0003
        paste0("OQ", v),                  # 752 -> OQ752
        paste0("OR", v),
        paste0("OQ", sprintf("%03d", as.integer(v))),
        paste0("OR", sprintf("%03d", as.integer(v)))
      )
    ))
    hit <- cand[cand %in% canonical]
    if (length(hit) == 0) NA_character_ else hit[1]
  }, character(1), USE.NAMES = FALSE)
}

# --- Source 1: identifier mapping (SRC0001) --------------------------------
id_code <- read.delim(src$id_code, stringsAsFactors = FALSE) %>%
  mutate(across(everything(), trimws)) %>%
  filter(nzchar(Code), nzchar(ID))
stopifnot(!any(duplicated(id_code$Code)), !any(duplicated(id_code$ID)))

# --- Cohort: the animals that actually have manual tracking data -----------
epm_raw   <- read.delim(src$epm,     stringsAsFactors = FALSE, check.names = FALSE)
epm_g_raw <- read.delim(src$epm_grp, stringsAsFactors = FALSE, check.names = FALSE)
nor_raw   <- suppressMessages(read_excel(src$nor))
socp_raw  <- suppressMessages(read_excel(src$socp))

cohort_codes <- sort(unique(c(trimws(epm_raw$Subject), trimws(nor_raw$Subject))))
cohort <- id_code %>% filter(Code %in% cohort_codes)
stopifnot(nrow(cohort) == length(cohort_codes))

# --- Source 2: curated group / sex assignment (RFID QC table) --------------
# NOTE: this table Group column holds the *phenotype* level (CON/RES/SUS),
# not the experimental condition. AnimalID_norm has lost zero-padding, so it
# is re-normalised before joining. Rows duplicated by that padding loss are
# collapsed; disagreement would show up as a pipe-joined value.
assign_raw <- read.csv(src$assign, stringsAsFactors = FALSE)
assign_tbl <- assign_raw %>%
  mutate(ID = normalise_id(AnimalID_norm, id_code$ID)) %>%
  filter(!is.na(ID)) %>%
  group_by(ID) %>%
  summarise(
    Batch_ref = paste(unique(na.omit(Batch_norm)),      collapse = "|"),
    Sex_ref   = paste(unique(na.omit(Sex)),             collapse = "|"),
    Pheno_ref = paste(unique(na.omit(Group)),           collapse = "|"),
    RefGroup  = paste(unique(na.omit(ReferenceGroup)),  collapse = "|"),
    n_rows    = n(),
    .groups   = "drop"
  )

# --- Source 3: expected phenotype (cross-scale QC) -------------------------
pheno_raw <- read.csv(src$phenotype, stringsAsFactors = FALSE)
pheno_tbl <- pheno_raw %>%
  mutate(ID = normalise_id(AnimalNum, id_code$ID)) %>%
  filter(!is.na(ID)) %>%
  group_by(ID) %>%
  summarise(Pheno_expected = paste(unique(na.omit(ExpectedGroup)), collapse = "|"),
            .groups = "drop")

# --- Sources 4-6: standalone SUS / CON animal lists ------------------------
#
# Analysis/.old/sus_animals.csv is deliberately NOT read: it is superseded.
# It is worth knowing it exists, though, because on the Batch-1 subset it
# matches sus_animals_batchCorrected.csv exactly (10 animals including 0001,
# OQ750 and OQ762), while the current sus_animals.csv drops those three. The
# batch-corrected call therefore reinstates the oldest list rather than
# inventing a new one. Recorded in README.md against the T2H7 conflict.
read_list <- function(p) {
  v <- trimws(readLines(p, warn = FALSE))
  normalise_id(v[nzchar(v)], id_code$ID)
}
sus_plain <- na.omit(read_list(src$sus))
sus_bc    <- na.omit(read_list(src$sus_bc))
con_list  <- na.omit(read_list(src$con))

# --- Condition (CON/SIS) as recorded in the assay scoring sheets -----------
cond_nor <- nor_raw %>%
  transmute(Code = trimws(Subject), Condition_NOR = trimws(Group)) %>%
  distinct()
cond_socp <- socp_raw %>%
  transmute(ID = normalise_id(ID, id_code$ID), Condition_SocP = trimws(Group)) %>%
  distinct() %>%
  filter(!is.na(ID))

# --- Assemble the canonical table ------------------------------------------
meta <- cohort %>%
  left_join(assign_tbl, by = "ID") %>%
  left_join(pheno_tbl,  by = "ID") %>%
  left_join(cond_nor,   by = "Code") %>%
  left_join(cond_socp,  by = "ID") %>%
  mutate(
    Batch = "B1",
    Sex   = ifelse(is.na(Sex_ref) | !nzchar(Sex_ref), NA_character_, Sex_ref),

    # Experimental condition, from the assay scoring sheets.
    Condition = coalesce(Condition_NOR, Condition_SocP),

    # Derived phenotype, primary definition: curated RFID QC table, falling
    # back to the standalone SUS/CON lists where that table has no row.
    Phenotype = case_when(
      !is.na(Pheno_ref) & nzchar(Pheno_ref) ~ Pheno_ref,
      ID %in% con_list                      ~ "CON",
      ID %in% sus_plain                     ~ "SUS",
      TRUE                                  ~ NA_character_
    ),
    Phenotype_source = case_when(
      !is.na(Pheno_ref) & nzchar(Pheno_ref) ~ "rfid_animal_group_sex_assignment_qc",
      ID %in% con_list                      ~ "Analysis/con_animals.csv",
      ID %in% sus_plain                     ~ "Analysis/sus_animals.csv",
      TRUE                                  ~ "unresolved"
    ),

    # Alternative phenotype definition: batch-corrected SUS call. Kept as its
    # own column because it disagrees with the primary definition for some
    # animals, and the choice changes downstream group comparisons.
    #
    # IMPORTANT ASYMMETRY: sus_animals_batchCorrected.csv lists only SUS
    # animals, so "RES" below is the *complement* within SIS, i.e. it assumes
    # every SIS animal absent from the list was scored and found resilient.
    # The primary definition makes no such assumption -- an animal with no row
    # in the curated table stays NA. Phenotype_bc_complement marks which RES
    # calls rest on that assumption so it can be excluded if unwarranted.
    Phenotype_batchCorrected = case_when(
      ID %in% sus_bc     ~ "SUS",
      Condition == "CON" ~ "CON",
      Condition == "SIS" ~ "RES",
      TRUE               ~ NA_character_
    ),
    Phenotype_bc_complement = Condition == "SIS" & !(ID %in% sus_bc),
    Phenotype_conflict = !is.na(Phenotype) & !is.na(Phenotype_batchCorrected) &
                         Phenotype != Phenotype_batchCorrected,

    # Condition implied by the curated phenotype, used as a consistency check.
    Condition_implied = case_when(
      Phenotype == "CON"             ~ "CON",
      Phenotype %in% c("RES", "SUS") ~ "SIS",
      TRUE                           ~ NA_character_
    ),
    Condition_conflict = !is.na(Condition) & !is.na(Condition_implied) &
                         Condition != Condition_implied,

    has_NOR  = Code %in% trimws(nor_raw$Subject),
    has_EPM  = Code %in% trimws(epm_raw$Subject),
    has_SocP = ID   %in% cond_socp$ID
  ) %>%
  select(Code, ID, Batch, Sex, Condition, Phenotype, Phenotype_source,
         Phenotype_batchCorrected, Phenotype_bc_complement, Phenotype_conflict,
         Condition_NOR, Condition_SocP, Condition_implied, Condition_conflict,
         Pheno_expected, RefGroup, has_NOR, has_EPM, has_SocP) %>%
  arrange(Code)

# --- Conflict / gap register -----------------------------------------------
conflicts <- bind_rows(
  meta %>% filter(Phenotype_conflict) %>%
    transmute(Code, ID, Domain = "phenotype",
              Description = sprintf(
                "Primary phenotype '%s' (%s) disagrees with batch-corrected call '%s' (Analysis/sus_animals_batchCorrected.csv).",
                Phenotype, Phenotype_source, Phenotype_batchCorrected),
              Status = "OPEN"),
  meta %>% filter(Condition_conflict) %>%
    transmute(Code, ID, Domain = "condition",
              Description = sprintf(
                "Assay-sheet condition '%s' disagrees with phenotype-implied condition '%s'.",
                Condition, Condition_implied),
              Status = "OPEN"),
  meta %>% filter(is.na(Phenotype)) %>%
    transmute(Code, ID, Domain = "phenotype",
              Description = "No RES/SUS phenotype found in any curated Exp9 source; left NA rather than imputed.",
              Status = "OPEN"),
  meta %>% filter(is.na(Sex)) %>%
    transmute(Code, ID, Domain = "sex",
              Description = "Not present in the curated RFID assignment table; Sex unknown from that source.",
              Status = "OPEN")
)

# --- Write metadata --------------------------------------------------------
write.table(meta, file.path(META, "animal_metadata.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write_xlsx(list(animal_metadata = meta, conflicts = conflicts),
           file.path(META, "animal_metadata.xlsx"))
write.csv(conflicts, file.path(META, "metadata_conflicts.csv"),
          row.names = FALSE, na = "")

# SLEAPanalyzer-compatible key files (tab separated, Code as the key), so the
# same metadata can be fed to the R pipeline via read_metadata_table().
write.table(meta %>% select(Code, ID), file.path(META, "animalIDCode.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
nor_novelloc <- nor_raw %>%
  transmute(Code = trimws(Subject), NovelLoc = trimws(novelloc)) %>%
  filter(!is.na(NovelLoc), nzchar(NovelLoc)) %>%
  distinct() %>%
  arrange(Code)
write.table(nor_novelloc, file.path(META, "novelLoc.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Provenance record -----------------------------------------------------
prov <- c(
  "Exp9 Batch-1 manual-tracking metadata - provenance",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("R:", R.version.string),
  "",
  "Sources (key | path | size bytes | last modified):",
  vapply(names(src), function(k) {
    fi <- file.info(src[[k]])
    sprintf("  %-9s %s | %s | %s", k, src[[k]], fi$size,
            format(fi$mtime, "%Y-%m-%d %H:%M:%S"))
  }, character(1)),
  "",
  sprintf("cohort: %d animals (union of EPM and NOR subjects)", nrow(meta)),
  sprintf("phenotype resolved: %d/%d", sum(!is.na(meta$Phenotype)), nrow(meta)),
  sprintf("phenotype conflicts (primary vs batch-corrected): %d", sum(meta$Phenotype_conflict)),
  sprintf("condition conflicts (assay sheet vs phenotype-implied): %d", sum(meta$Condition_conflict)),
  "",
  "Note: the curated RFID table Group column encodes phenotype (CON/RES/SUS),",
  "not experimental condition. Condition here comes from the assay scoring sheets.",
  "Original NOR.xlsx / SocP.xlsx / *.tsv files were read only, never written."
)
writeLines(prov, file.path(META, "metadata_provenance.txt"))

cat("\n=== canonical metadata ===\n")
print(as.data.frame(meta %>% select(Code, ID, Sex, Condition, Phenotype,
                                    Phenotype_batchCorrected, Phenotype_bc_complement, Phenotype_conflict,
                                    has_NOR, has_EPM, has_SocP)))
cat("\n=== conflicts / gaps ===\n")
print(as.data.frame(conflicts))
cat("\nwrote metadata to:", META, "\n")
