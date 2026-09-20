# ============================================================================
# 04_assemble_sleap.R
#
# Reads the SLEAPanalyzer v2 output for EPM, NOR and SocP (S1, S2), joins the
# canonical metadata, and writes one per-animal SLEAP table to match
# enriched/analysis_ready_wide.tsv on `Code`.
#
# Conventions that matter when reading this table:
#
#   NOR sides are MIRRORED relative to the manual BORIS scoring. SLEAP's objL
#   is image-left in all 20 files (objL_x < objR_x), and SLEAP's tl is left of
#   tr, so the SLEAP export is internally consistent; it is the manual
#   intLeft / intRight columns that correspond to the opposite image side.
#   Evidence: contactLeft ~ intRight_BORIS r = +0.96 while
#   contactLeft ~ intLeft_BORIS r = +0.29.
#
#   The pipeline's novel/familiar mapping (metadata "R" -> novel is LEFT,
#   Behavioral_Metrics_Phase1.R:212) compensates for exactly that mirror, so
#   contactNov / contactFam DO align with nov_BORIS / fam_BORIS. Use those,
#   not the side columns.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(writexl)
})

PROJ   <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
OUT    <- file.path(PROJ, "sleap_output")
ENRICH <- file.path(PROJ, "enriched")

meta <- read.delim(file.path(PROJ, "metadata/animal_metadata.tsv"),
                   stringsAsFactors = FALSE, na.strings = "")

code_of <- function(x) substr(basename(x), 1, 4)

# ---------------------------------------------------------------- EPM -------
epm <- read.csv(file.path(OUT, "EPM/Report.csv"), stringsAsFactors = FALSE) %>%
  mutate(Code = code_of(file)) %>%
  select(
    Code,
    sleap_EPM_open_s    = bodycentre.open.total.time,
    sleap_EPM_closed_s  = bodycentre.closed.total.time,
    sleap_EPM_center_s  = bodycentre.center.total.time,
    sleap_EPM_total_s   = bodycentre.total.time,
    sleap_EPM_open_entries   = bodycentre.open.entries,
    sleap_EPM_closed_entries = bodycentre.closed.entries,
    sleap_EPM_center_entries = bodycentre.center.entries,
    sleap_EPM_nosedip_n = nose.dip,
    sleap_EPM_distance_cm = bodycentre.raw.distance,
    sleap_EPM_speed_raw   = bodycentre.raw.speed,
    sleap_EPM_speed_moving = bodycentre.speed.moving,
    sleap_EPM_time_moving  = bodycentre.time.moving,
    sleap_EPM_time_stationary = bodycentre.time.stationary
  ) %>%
  mutate(
    sleap_EPM_open_frac   = sleap_EPM_open_s / sleap_EPM_total_s,
    sleap_EPM_closed_frac = sleap_EPM_closed_s / sleap_EPM_total_s
  )

# ---------------------------------------------------------------- NOR -------
nor <- read.csv(file.path(OUT, "NOR/combined_output.csv"), stringsAsFactors = FALSE) %>%
  mutate(Code = code_of(file)) %>%
  select(
    Code,
    sleap_NOR_nov_s   = contactNov,
    sleap_NOR_fam_s   = contactFam,
    sleap_NOR_left_s  = contactLeft,
    sleap_NOR_right_s = contactRight,
    sleap_NOR_proxNov_s = proxNov,
    sleap_NOR_proxFam_s = proxFam,
    sleap_NOR_latency_s = latency,
    sleap_NOR_freqL = frequencyL,
    sleap_NOR_freqR = frequencyR,
    sleap_NOR_distance_cm = distance,
    sleap_NOR_speedMoving = speedMoving,
    sleap_NOR_rear_n = frequencyRear_experimental,
    sleap_NOR_novelLoc = novelLoc,
    sleap_NOR_qcPass = qcPass
  ) %>%
  mutate(
    sleap_NOR_total_s = sleap_NOR_nov_s + sleap_NOR_fam_s,
    sleap_NOR_D2 = (sleap_NOR_nov_s - sleap_NOR_fam_s) / sleap_NOR_total_s
  )

# --------------------------------------------------------------- SocP -------
socp <- bind_rows(lapply(c("S1", "S2"), function(ph) {
  read.csv(file.path(OUT, "SocP", ph, "combined_output.csv"), stringsAsFactors = FALSE) %>%
    mutate(Code = code_of(file), phase = ph)
})) %>%
  select(Code, phase,
         novel = contactNovel, familiar = contactFamiliar,
         proxNovel, proxFamiliar,
         left = contactLeft, right = contactRight,
         entriesLeft, entriesRight) %>%
  mutate(
    total = novel + familiar,
    pref_index = (novel - familiar) / total
  ) %>%
  pivot_wider(
    names_from  = phase,
    values_from = c(novel, familiar, proxNovel, proxFamiliar, left, right,
                    entriesLeft, entriesRight, total, pref_index),
    names_glue  = "sleap_SocP_{phase}_{.value}"
  )

# ------------------------------------------------------------- assemble -----
sleap <- meta %>%
  select(Code, ID, Batch, Sex, Condition, Phenotype, Phenotype_batchCorrected,
         Phenotype_bc_complement, Phenotype_conflict) %>%
  left_join(epm,  by = "Code") %>%
  left_join(nor,  by = "Code") %>%
  left_join(socp, by = "Code") %>%
  arrange(Code)

stopifnot(nrow(sleap) == nrow(meta), !any(duplicated(sleap$Code)))

write.table(sleap, file.path(ENRICH, "sleap_wide.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write_xlsx(sleap, file.path(ENRICH, "sleap_wide.xlsx"))

# ------------------------------------------------------------- report -------
cat(sprintf("sleap_wide: %d animals x %d columns\n", nrow(sleap), ncol(sleap)))
num <- names(sleap)[vapply(sleap, is.numeric, logical(1))]
miss <- vapply(sleap[num], function(x) sum(is.na(x)), integer(1))
cat(sprintf("numeric metrics: %d | columns with NA: %d\n", length(num), sum(miss > 0)))
if (any(miss > 0)) print(miss[miss > 0])

cat("\n=== magnitude comparison: is each SLEAP/BORIS pair the same measure? ===\n")
boris <- read.delim(file.path(ENRICH, "analysis_ready_wide.tsv"),
                    stringsAsFactors = FALSE, na.strings = "")
cmp <- function(label, s, b) {
  cat(sprintf("  %-26s SLEAP %7.1f [%6.1f-%7.1f] | BORIS %7.1f [%6.1f-%7.1f] | ratio %.2f\n",
              label, median(s), min(s), max(s), median(b), min(b), max(b),
              median(s) / median(b)))
}
j <- inner_join(sleap, boris, by = "Code")
cmp("EPM open time (s)",     j$sleap_EPM_open_s,   j$EPM_OpenTime_dur)
cmp("EPM closed time (s)",   j$sleap_EPM_closed_s, j$EPM_ClosedTime_dur)
cmp("EPM center time (s)",   j$sleap_EPM_center_s, j$EPM_CenterTime_dur)
cmp("EPM nose dips (count)", j$sleap_EPM_nosedip_n, j$EPM_ND_n)
cmp("NOR novel (s)",         j$sleap_NOR_nov_s,    j$NOR_nov)
cmp("NOR familiar (s)",      j$sleap_NOR_fam_s,    j$NOR_fam)
cmp("SocP S1 novel (s)",     j$sleap_SocP_S1_novel, j$SocP_S1_novel)
cmp("SocP S1 familiar (s)",  j$sleap_SocP_S1_familiar, j$SocP_S1_familiar)
cmp("SocP S2 novel (s)",     j$sleap_SocP_S2_novel, j$SocP_S2_novel)
cmp("SocP S2 familiar (s)",  j$sleap_SocP_S2_familiar, j$SocP_S2_familiar)

cat("\n=== NOR novel/familiar assignment check ===\n")
cat(sprintf("  contactNov ~ nov_BORIS   r = %+.3f   (expected to be the match)\n",
            cor(j$sleap_NOR_nov_s, j$NOR_nov)))
cat(sprintf("  contactNov ~ fam_BORIS   r = %+.3f\n", cor(j$sleap_NOR_nov_s, j$NOR_fam)))
cat(sprintf("  D2:  sleap ~ boris       r = %+.3f\n", cor(j$sleap_NOR_D2, j$NOR_D2)))

cat("\n=== SocP novel/familiar assignment check (is it mirrored too?) ===\n")
for (ph in c("S1", "S2")) {
  sn <- j[[paste0("sleap_SocP_", ph, "_novel")]]
  bn <- j[[paste0("SocP_", ph, "_novel")]]
  bf <- j[[paste0("SocP_", ph, "_familiar")]]
  sp <- j[[paste0("sleap_SocP_", ph, "_pref_index")]]
  bp <- j[[paste0("SocP_", ph, "_pref_index")]]
  cat(sprintf("  %s: novel~novel r = %+.3f | novel~familiar r = %+.3f | pref_index r = %+.3f\n",
              ph, cor(sn, bn), cor(sn, bf), cor(sp, bp)))
}
cat("\nwrote:", file.path(ENRICH, "sleap_wide.tsv"), "\n")
