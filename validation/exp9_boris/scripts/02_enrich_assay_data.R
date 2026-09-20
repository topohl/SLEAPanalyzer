# ============================================================================
# 02_enrich_assay_data.R
#
# Joins the canonical Exp9 metadata (built by 01_build_metadata.R) onto the
# three manual-scoring files, and assembles one analysis-ready per-animal
# table for the SLEAP-vs-BORIS correlation work.
#
# Originals are read-only. Everything is written to enriched/.
#
# Conventions worth knowing before interpreting the output:
#
#   NOR novel/familiar mapping. The manual sheet uses the DIRECT convention:
#   novelloc == "L" means nov_BORIS == intLeft_BORIS. SLEAPanalyzer's
#   documented (inherited) convention is the INVERSE -- docs/assay_definitions.md
#   states that metadata "R" maps the LEFT side to novel. Feeding the same
#   novelLoc table to both therefore swaps novel and familiar between the
#   manual and automated data. See metadata/novelLoc.txt and the note below.
#
#   D2 = (nov - fam) / (nov + fam), verified against D2_BORIS for all 20.
#
#   BatchZscoreD2_BORIS is not a batch z-score: it equals
#   (D2 - mean(D2_CON)) / sd_pop(D2_CON), i.e. referenced to the 4 CON
#   animals using the population SD (denominator n, not n-1).
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(writexl)
  library(stringr)
})

PROJ   <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
META   <- file.path(PROJ, "metadata")
ENRICH <- file.path(PROJ, "enriched")
dir.create(ENRICH, showWarnings = FALSE, recursive = TRUE)

meta <- read.delim(file.path(META, "animal_metadata.tsv"),
                   stringsAsFactors = FALSE, na.strings = "")
meta_join <- meta %>%
  select(Code, ID, Batch, Sex, Condition, Phenotype, Phenotype_source,
         Phenotype_batchCorrected, Phenotype_bc_complement, Phenotype_conflict)

write_both <- function(df, stem) {
  write.table(df, file.path(ENRICH, paste0(stem, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  write_xlsx(df, file.path(ENRICH, paste0(stem, ".xlsx")))
  cat(sprintf("  %-28s %3d x %2d\n", stem, nrow(df), ncol(df)))
}

cat("writing enriched files:\n")

# ============================ NOR ==========================================
nor_raw <- suppressMessages(read_excel(file.path(PROJ, "NOR.xlsx")))

nor <- nor_raw %>%
  rename(ID_sheet = ID, Code = Subject, Group_sheet = Group) %>%
  mutate(Code = trimws(Code)) %>%
  left_join(meta_join, by = "Code") %>%
  mutate(
    # Recomputed discrimination index, and the side-mapping audit.
    NOR_D2_recomputed = (nov_BORIS - fam_BORIS) / (nov_BORIS + fam_BORIS),
    NOR_D2_ok = abs(NOR_D2_recomputed - D2_BORIS) < 1e-6,

    # Does the novel/familiar pair reconcile with the left/right pair? If not,
    # the two column pairs came from different scoring passes.
    NOR_sides_reconcile =
      (abs(nov_BORIS - intLeft_BORIS)  < 1e-9 & abs(fam_BORIS - intRight_BORIS) < 1e-9) |
      (abs(nov_BORIS - intRight_BORIS) < 1e-9 & abs(fam_BORIS - intLeft_BORIS)  < 1e-9),

    # Which side the sheet actually treated as novel, read back from the data.
    NOR_novel_side_observed = case_when(
      abs(nov_BORIS - intLeft_BORIS)  < 1e-9 ~ "L",
      abs(nov_BORIS - intRight_BORIS) < 1e-9 ~ "R",
      TRUE                                   ~ NA_character_
    ),
    NOR_mapping_convention = ifelse(
      !is.na(NOR_novel_side_observed) & NOR_novel_side_observed == trimws(novelloc),
      "direct", "inverted_or_unresolved"
    ),
    NOR_total_interaction = nov_BORIS + fam_BORIS
  )

# CON-referenced z-score, reproducing the sheet's BatchZscoreD2_BORIS, plus a
# conventional (n-1) version. Both rest on only 4 CON animals, so they are
# extremely unstable -- reported for traceability, not as primary outcomes.
con_d2 <- nor$D2_BORIS[nor$Condition == "CON"]
sd_pop <- function(x) sqrt(mean((x - mean(x))^2))
nor <- nor %>%
  mutate(
    NOR_D2_zCON_pop    = (D2_BORIS - mean(con_d2)) / sd_pop(con_d2),
    NOR_D2_zCON_sample = (D2_BORIS - mean(con_d2)) / sd(con_d2),
    NOR_zCON_n_ref     = length(con_d2)
  ) %>%
  relocate(Code, ID, ID_sheet, Batch, Sex, Condition, Phenotype)

write_both(nor, "NOR_enriched")

# ============================ SocP =========================================
socp_raw <- suppressMessages(read_excel(file.path(PROJ, "SocP.xlsx")))

# SocP.xlsx carries only a numeric ID (Excel dropped the OQ prefix and the
# leading zeros), so it is joined through the canonical ID reconstructed the
# same way 01_build_metadata.R does it.
socp_long <- socp_raw %>%
  rename(ID_sheet = ID, Group_sheet = Group) %>%
  mutate(
    ID = vapply(ID_sheet, function(v) {
      v <- trimws(as.character(v))
      cand <- c(v, sprintf("%04d", suppressWarnings(as.integer(v))), paste0("OQ", v))
      hit <- cand[cand %in% meta$ID]
      if (length(hit) == 0) NA_character_ else hit[1]
    }, character(1), USE.NAMES = FALSE)
  ) %>%
  left_join(meta_join, by = "ID") %>%
  relocate(Code, ID, ID_sheet, Batch, Sex, Condition, Phenotype)
stopifnot(!any(is.na(socp_long$Code)))

write_both(socp_long, "SocP_enriched_long")

# Per-animal wide form with a preference index per phase.
socp_wide <- socp_long %>%
  select(Code, Side,
         S1 = `S1Total duration (s)`,
         S2 = `S2Total duration (s)`) %>%
  pivot_longer(c(S1, S2), names_to = "Phase", values_to = "duration") %>%
  pivot_wider(names_from = Side, values_from = duration) %>%
  mutate(
    total = novel + familiar,
    pref_index = (novel - familiar) / total
  ) %>%
  pivot_wider(
    names_from  = Phase,
    values_from = c(novel, familiar, total, pref_index),
    names_glue  = "SocP_{Phase}_{.value}"
  )

write_both(socp_wide, "SocP_wide")

# ============================ EPM ==========================================
epm_raw <- read.delim(file.path(PROJ, "SISB1_EPMdata.tsv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
epm_g_raw <- read.delim(file.path(PROJ, "SISB1_EPMdataGrouped.tsv"),
                        stringsAsFactors = FALSE, check.names = FALSE)

# BORIS writes NA, not 0, into the duration and percentage columns when a
# behaviour was never observed (Total number of occurences == 0). Those are
# structural zeros: the animal was observed and did not perform the behaviour.
# Left as NA they would silently drop animals from every correlation that
# touches the column, so they are recoded to 0 and the recode is counted.
epm_structural_zeros <- epm_raw %>%
  filter(`Total number of occurences` == 0) %>%
  count(Behavior, name = "n_animals_zero")

epm_raw <- epm_raw %>%
  mutate(
    epm_structural_zero = `Total number of occurences` == 0,
    `Total duration (s)` = ifelse(epm_structural_zero, 0, `Total duration (s)`),
    `% of total length`  = ifelse(epm_structural_zero, 0, `% of total length`)
  )

epm_long <- epm_raw %>%
  mutate(Code = trimws(Subject)) %>%
  left_join(meta_join, by = "Code") %>%
  relocate(Code, ID, Batch, Sex, Condition, Phenotype)
stopifnot(!any(is.na(epm_long$ID)))

epm_g_long <- epm_g_raw %>%
  mutate(`Total duration (s)` = ifelse(`Total number` == 0, 0, `Total duration (s)`)) %>%
  mutate(Code = trimws(Subject)) %>%
  left_join(meta_join, by = "Code") %>%
  relocate(Code, ID, Batch, Sex, Condition, Phenotype)

write_both(epm_long, "EPM_enriched_long")
write_both(epm_g_long, "EPM_grouped_enriched_long")

# Per-animal wide form: duration and count for every scored behaviour.
epm_beh <- epm_long %>%
  select(Code, Behavior,
         dur = `Total duration (s)`,
         n   = `Total number of occurences`,
         pct = `% of total length`) %>%
  pivot_wider(names_from = Behavior, values_from = c(dur, n, pct),
              names_glue = "EPM_{Behavior}_{.value}")

epm_zone <- epm_g_long %>%
  select(Code, Category,
         dur = `Total duration (s)`,
         n   = `Total number`) %>%
  pivot_wider(names_from = Category, values_from = c(dur, n),
              names_glue = "EPM_{Category}_{.value}")

epm_wide <- epm_zone %>%
  left_join(epm_beh, by = "Code") %>%
  mutate(
    EPM_arena_time   = EPM_OpenTime_dur + EPM_CenterTime_dur + EPM_ClosedTime_dur,
    EPM_open_frac    = EPM_OpenTime_dur / EPM_arena_time,
    EPM_closed_frac  = EPM_ClosedTime_dur / EPM_arena_time,
    EPM_total_entries = EPM_OpenTime_n + EPM_CenterTime_n + EPM_ClosedTime_n,
    # Risk assessment and affective-state composites, summed across zones.
    EPM_SAP_dur   = EPM_OpenSAP_dur + EPM_CenterSAP_dur + EPM_ClosedSAP_dur,
    EPM_SAP_n     = EPM_OpenSAP_n   + EPM_CenterSAP_n   + EPM_ClosedSAP_n,
    EPM_ND_dur    = EPM_OpenND_dur  + EPM_CenterND_dur,
    EPM_ND_n      = EPM_OpenND_n    + EPM_CenterND_n,
    EPM_groom_dur = EPM_OpenGroom_dur + EPM_ClosedGroom_dur,
    EPM_rear_dur  = EPM_CenterRear_dur + EPM_ClosedRear_dur,
    EPM_explor_dur = EPM_OpenExplor_dur + EPM_CenterExplor_dur + EPM_ClosedExplor_dur
  )

write_both(epm_wide, "EPM_wide")

# ==================== Analysis-ready master table ==========================
nor_wide <- nor %>%
  select(Code,
         NOR_intLeft  = intLeft_BORIS,
         NOR_intRight = intRight_BORIS,
         NOR_nov      = nov_BORIS,
         NOR_fam      = fam_BORIS,
         NOR_D2       = D2_BORIS,
         NOR_total_interaction,
         NOR_D2_zCON_pop,
         NOR_novelloc = novelloc,
         NOR_sides_reconcile,
         NOR_mapping_convention)

master <- meta %>%
  select(Code, ID, Batch, Sex, Condition, Phenotype, Phenotype_source,
         Phenotype_batchCorrected, Phenotype_bc_complement, Phenotype_conflict) %>%
  left_join(nor_wide,  by = "Code") %>%
  left_join(socp_wide, by = "Code") %>%
  left_join(epm_wide,  by = "Code") %>%
  arrange(Code)

stopifnot(nrow(master) == nrow(meta), !any(duplicated(master$Code)))

write_both(master, "analysis_ready_wide")

# ============================ Report =======================================
cat("\n=== master table: metric coverage ===\n")
num_cols <- names(master)[vapply(master, is.numeric, logical(1))]
cat(sprintf("%d animals, %d numeric metrics, %d columns total\n",
            nrow(master), length(num_cols), ncol(master)))
miss <- vapply(master[num_cols], function(x) sum(is.na(x)), integer(1))
if (any(miss > 0)) {
  cat("numeric columns with missing values:\n"); print(miss[miss > 0])
} else {
  cat("no missing values in any numeric metric\n")
}

cat("\n=== EPM structural zeros recoded (count == 0 -> duration 0) ===\n")
print(as.data.frame(epm_structural_zeros))

# A metric that is zero for nearly every animal carries almost no information,
# and a Pearson correlation on it is driven entirely by the one or two animals
# that differ. Flag those rather than silently correlating them.
cat("\n=== low-variance metrics (unusable or fragile for correlation) ===\n")
variability <- data.frame(
  metric   = num_cols,
  n_unique = vapply(master[num_cols], function(x) length(unique(na.omit(x))), integer(1)),
  n_nonzero = vapply(master[num_cols], function(x) sum(na.omit(x) != 0), integer(1)),
  sd       = vapply(master[num_cols], function(x) sd(x, na.rm = TRUE), numeric(1)),
  row.names = NULL
)
fragile <- variability %>% filter(n_unique <= 3 | n_nonzero <= 3) %>% arrange(n_nonzero)
if (nrow(fragile) > 0) print(fragile) else cat("none\n")
write.csv(variability, file.path(ENRICH, "metric_variability.csv"), row.names = FALSE)

cat("\n=== NOR side-mapping audit ===\n")
print(nor %>% count(novelloc, NOR_novel_side_observed, NOR_mapping_convention))
cat("rows where nov/fam do not reconcile with left/right:\n")
print(as.data.frame(nor %>% filter(!NOR_sides_reconcile) %>%
  select(Code, ID, novelloc, intLeft_BORIS, intRight_BORIS, nov_BORIS, fam_BORIS, D2_BORIS)))

cat("\n=== D2 reproduction check ===\n")
cat("D2_BORIS reproduced by (nov-fam)/(nov+fam):", all(nor$NOR_D2_ok), "\n")
cat("BatchZscoreD2_BORIS reproduced by CON-referenced population z:",
    isTRUE(all.equal(nor$NOR_D2_zCON_pop, nor$BatchZscoreD2_BORIS)),
    sprintf("(n_CON = %d)", unique(nor$NOR_zCON_n_ref)), "\n")

cat("\n=== group sizes ===\n")
cat("Condition:\n");                print(table(master$Condition, useNA = "ifany"))
cat("Phenotype (primary):\n");      print(table(master$Phenotype, useNA = "ifany"))
cat("Phenotype (batch-corrected):\n"); print(table(master$Phenotype_batchCorrected, useNA = "ifany"))

cat("\nwrote enriched files to:", ENRICH, "\n")
