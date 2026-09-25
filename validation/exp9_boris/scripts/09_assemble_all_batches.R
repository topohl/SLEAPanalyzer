# ============================================================================
# 09_assemble_all_batches.R
#
# Builds the full-cohort metadata table and assembles the B1-B6 SLEAP output
# for EPM, NOR, SocP and OFT into one per-animal table.
#
# 01_build_metadata.R deliberately restricts itself to the 20 Batch-1 animals
# that have manual BORIS scoring. This script applies the same sources and the
# same rules to every animal that appears in any SLEAP output, so the cohort
# is defined by the tracking data rather than by the manual scoring.
#
# Phenotype uses only the corrected canonical lists in Analysis/: controls are
# named by con_animals.csv, susceptible animals by sus_animals.txt, and RES is
# the complement among SIS animals. Archived alternative rosters are excluded.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(writexl)
})

EXP9 <- Sys.getenv(
  "EXP9_ROOT",
  unset = "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress"
)
LEGACY_PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
RUN_ROOT <- Sys.getenv("EXP9_SLEAP_RUN_ROOT", unset = LEGACY_PROJ)
STAGE <- Sys.getenv(
  "EXP9_SLEAP_STAGE_ROOT",
  unset = "C:/Users/topohl/Documents/exp9_sleap_staging"
)
OUT <- Sys.getenv(
  "EXP9_SLEAP_ASSAY_OUTPUT_ROOT",
  unset = file.path(LEGACY_PROJ, "sleap_output_all")
)
ENRICH <- Sys.getenv("EXP9_SLEAP_DATA_DIR", unset = file.path(RUN_ROOT, "enriched"))
META <- Sys.getenv("EXP9_SLEAP_METADATA_DIR", unset = file.path(RUN_ROOT, "metadata"))
BATCHES <- paste0("B", 1:6)
dir.create(ENRICH, showWarnings = FALSE, recursive = TRUE)
dir.create(META, showWarnings = FALSE, recursive = TRUE)

# --- Metadata for the whole cohort -----------------------------------------
id_code <- read.delim(file.path(EXP9, "Planning/animalIDCode.txt"),
                      stringsAsFactors = FALSE) %>%
  mutate(across(everything(), trimws)) %>%
  filter(nzchar(Code), nzchar(ID))

# Exp9 writes the same animal several ways across files, and crucially the
# RFID table's AnimalID_norm STRIPS LEADING ZEROS (00318 -> 318) while
# animalIDCode.txt keeps them. Matching on the literal strings silently lost
# 21 animals. Reduce both sides to one canonical form instead: strip leading
# zeros from purely numeric identifiers, upper-case everything else. No
# numeric conversion, so long identifiers cannot overflow to NA.
canon <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(grepl("^[0-9]+$", x), sub("^0+(?=[0-9])", "", x, perl = TRUE), x)
}

id_code$key <- canon(id_code$ID)

# The RFID table also carries two identifier systems per row: AnimalID_raw is
# the one animalIDCode uses for some batches, AnimalID_norm for others. Accept
# a match on either. Read the Stage 01 foundation copy, which is identical to
# the numbered 03_derived_metrics original that is being archived.
asg_raw <- read.csv(file.path(EXP9, "Analysis/Behavior/RFID/analysis_ready/foundations/behavior_metrics/qc",
                              "animal_group_sex_assignment_qc.csv"), stringsAsFactors = FALSE)
asg <- bind_rows(
    asg_raw %>% mutate(key = canon(AnimalID_norm)),
    asg_raw %>% mutate(key = canon(AnimalID_raw))
  ) %>%
  filter(key %in% id_code$key) %>%
  group_by(key) %>%
  summarise(Batch_ref  = paste(unique(na.omit(Batch_norm)), collapse = "|"),
            Sex_rfid   = paste(unique(na.omit(Sex)), collapse = "|"),
            Pheno_rfid = paste(unique(na.omit(Group)), collapse = "|"),
            .groups = "drop")

read_list <- function(p) {
  v <- trimws(readLines(p, warn = FALSE))
  canon(v[nzchar(v)])
}
# sus_animals.txt is the authoritative, corrected SUS list: it names every
# susceptible animal, and RES is the complement within SIS. The historical
# batch-corrected alternative was archived on 2026-09-20 as non-canonical and
# must not be silently reintroduced here.
sus_plain <- read_list(file.path(EXP9, "Analysis/sus_animals.txt"))
con_list  <- read_list(file.path(EXP9, "Analysis/con_animals.csv"))

# Sex is a property of the batch: B1, B2 and B5 are male, B3, B4 and B6 female.
# No batch is mixed, so sex is fully determined even for animals the RFID table
# never reached -- and equally, a sex effect is inseparable from a batch effect.
SEX_BY_BATCH <- c(B1 = "Male", B2 = "Male", B3 = "Female",
                  B4 = "Female", B5 = "Male", B6 = "Female")

# --- Assay readers ----------------------------------------------------------
per_batch <- function(path_fun, tag) {
  bind_rows(lapply(BATCHES, function(b) {
    p <- path_fun(b)
    if (!file.exists(p)) return(NULL)
    read.csv(p, stringsAsFactors = FALSE) %>% mutate(batch = b)
  }))
}

epm <- per_batch(function(b) file.path(OUT, "EPM", b, "Report.csv")) %>%
  mutate(Code = substr(basename(file), 1, 4)) %>%
  select(Code, batch,
         sleap_EPM_open_s = bodycentre.open.total.time,
         sleap_EPM_closed_s = bodycentre.closed.total.time,
         sleap_EPM_center_s = bodycentre.center.total.time,
         sleap_EPM_total_s = bodycentre.total.time,
         sleap_EPM_open_entries = bodycentre.open.entries,
         sleap_EPM_nosedip_n = nose.dip,
         sleap_EPM_distance_cm = bodycentre.raw.distance,
         sleap_EPM_speed_moving = bodycentre.speed.moving,
         sleap_EPM_time_moving = bodycentre.time.moving) %>%
  mutate(sleap_EPM_open_frac = sleap_EPM_open_s / sleap_EPM_total_s)

nor <- per_batch(function(b) file.path(OUT, "NOR", b, "combined_output.csv")) %>%
  select(Code, batch,
         sleap_NOR_nov_s = contactNov, sleap_NOR_fam_s = contactFam,
         sleap_NOR_left_s = contactLeft, sleap_NOR_right_s = contactRight,
         sleap_NOR_latency_s = latency, sleap_NOR_distance_cm = distance,
         sleap_NOR_novelLoc = novelLoc, sleap_NOR_qcPass = qcPass) %>%
  mutate(sleap_NOR_total_s = sleap_NOR_nov_s + sleap_NOR_fam_s,
         sleap_NOR_D2 = (sleap_NOR_nov_s - sleap_NOR_fam_s) / sleap_NOR_total_s)

socp <- bind_rows(lapply(BATCHES, function(b) {
  bind_rows(lapply(c("S1", "S2"), function(ph) {
    p <- file.path(OUT, "SocP", b, ph, "combined_output.csv")
    if (!file.exists(p)) return(NULL)
    read.csv(p, stringsAsFactors = FALSE) %>% mutate(batch = b, ph = ph)
  }))
})) %>%
  select(Code, batch, ph, novel = contactNovel, familiar = contactFamiliar,
         left = contactLeft, right = contactRight) %>%
  mutate(total = novel + familiar,
         pref_index = (novel - familiar) / total) %>%
  pivot_wider(names_from = ph,
              values_from = c(novel, familiar, left, right, total, pref_index),
              names_glue = "sleap_SocP_{ph}_{.value}")

oft_files <- unlist(lapply(BATCHES, function(b)
  list.files(file.path(STAGE, "OFT", b, "output_v1.2.0", "tables"),
             pattern = "_summary[.]csv$", full.names = TRUE)))
# ID arrives as character in some per-animal summaries ("0001") and integer in
# others ("750"), so bind_rows() refuses to combine them. Coerce the identifier
# columns before binding rather than letting the type depend on which animals
# happen to be in a batch.
oft <- if (length(oft_files) == 0) NULL else
  bind_rows(lapply(oft_files, function(p) {
    read.csv(p, stringsAsFactors = FALSE) %>%
      mutate(across(any_of(c("ID", "Code", "Batch", "file")), as.character))
  }))

cat(sprintf("assay rows read -- EPM %d | NOR %d | SocP %d | OFT %s\n",
            nrow(epm), nrow(nor), nrow(socp),
            if (is.null(oft)) "0 (not available)" else nrow(oft)))

# --- Cohort and metadata ----------------------------------------------------
cohort_codes <- sort(unique(c(epm$Code, nor$Code, socp$Code,
                              if (!is.null(oft)) oft$Code)))
cohort_codes <- cohort_codes[!is.na(cohort_codes)]

batch_of <- bind_rows(
  epm %>% select(Code, batch), nor %>% select(Code, batch),
  socp %>% select(Code, batch)) %>%
  distinct(Code, batch) %>% group_by(Code) %>%
  summarise(Batch = paste(sort(unique(batch)), collapse = "|"), .groups = "drop")

meta <- id_code %>%
  filter(Code %in% cohort_codes) %>%
  left_join(asg, by = "key") %>%
  left_join(batch_of, by = "Code") %>%
  mutate(
    # Condition is a roster fact: con_animals.csv names the controls (4 per
    # batch, 24 in total) and everything else in the cohort is stressed.
    Condition = ifelse(key %in% con_list, "CON", "SIS"),

    # Phenotype follows the lab's own convention: sus_animals.txt lists every
    # susceptible animal and RES is the complement within SIS. There is no
    # "unresolved" category -- absence from the list IS the resilient call.
    Phenotype = case_when(Condition == "CON" ~ "CON",
                          key %in% sus_plain ~ "SUS",
                          TRUE               ~ "RES"),
    Phenotype_source = "Analysis/{con_animals.csv, sus_animals.txt}, RES by complement",

    # The RFID QC table is an independent third source. It is NOT used to
    # assign anything -- it cross-checks the rosters, and disagreements are
    # reported rather than silently preferred either way.
    Pheno_rfid = ifelse(is.na(Pheno_rfid) | !nzchar(Pheno_rfid), NA_character_, Pheno_rfid),
    Phenotype_rfid_agrees = ifelse(is.na(Pheno_rfid), NA, Pheno_rfid == Phenotype),

    # Sex is determined by batch; the RFID value is kept only to verify that.
    Sex_rfid = ifelse(is.na(Sex_rfid) | !nzchar(Sex_rfid), NA_character_, Sex_rfid),
    Sex = unname(SEX_BY_BATCH[Batch]),
    Sex_rfid_agrees = ifelse(is.na(Sex_rfid), NA, Sex_rfid == Sex)) %>%
  select(Code, ID, Batch, Batch_ref, Sex, Sex_rfid, Sex_rfid_agrees,
         Condition, Phenotype, Phenotype_source,
         Pheno_rfid, Phenotype_rfid_agrees)

oft_wide <- if (is.null(oft)) NULL else
  oft %>%
    select(Code,
           sleap_OFT_distance_cm = distance_cm,
           sleap_OFT_speed_moving = moving_speed_cm_s,
           sleap_OFT_percent_moving = percent_moving,
           sleap_OFT_center_s = center_time_s,
           sleap_OFT_center_pct = center_time_percent,
           sleap_OFT_center_entries = center_entries,
           sleap_OFT_center_latency_s = center_latency_s,
           sleap_OFT_periphery_pct = periphery_time_percent,
           sleap_OFT_corner_pct = corner_time_percent,
           sleap_OFT_immobility_s = immobility_time_s,
           sleap_OFT_wall_distance_cm = mean_wall_distance_cm) %>%
    distinct(Code, .keep_all = TRUE)

master <- meta %>%
  left_join(epm  %>% select(-batch), by = "Code") %>%
  left_join(nor  %>% select(-batch), by = "Code") %>%
  left_join(socp %>% select(-batch), by = "Code") %>%
  {if (is.null(oft_wide)) . else left_join(., oft_wide, by = "Code")} %>%
  arrange(Batch, Code)
stopifnot(!any(duplicated(master$Code)))

write.table(master, file.path(ENRICH, "sleap_all_batches_wide.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write_xlsx(master, file.path(ENRICH, "sleap_all_batches_wide.xlsx"))
write.table(meta, file.path(META, "animal_metadata_all.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# --- Report -----------------------------------------------------------------
cat(sprintf("\ncohort: %d animals x %d columns\n", nrow(master), ncol(master)))

cat("\n=== assay coverage by batch ===\n")
print(master %>%
  mutate(has_EPM = !is.na(sleap_EPM_open_s), has_NOR = !is.na(sleap_NOR_nov_s),
         has_SocP = !is.na(sleap_SocP_S1_novel), has_OFT = !is.na(sleap_OFT_center_s)) %>%
  group_by(Batch) %>%
  # Distinct names from the summarised ones: summarise() evaluates
  # sequentially, so reusing EPM would make all_three read the scalar sum.
  summarise(animals = n(), EPM = sum(has_EPM), NOR = sum(has_NOR),
            SocP = sum(has_SocP), OFT = sum(has_OFT),
            all_four = sum(has_EPM & has_NOR & has_SocP & has_OFT), .groups = "drop") %>%
  as.data.frame(), row.names = FALSE)

cat("\n=== metadata coverage ===\n")
cat(sprintf("  Sex known            : %d / %d\n", sum(!is.na(master$Sex)), nrow(master)))
cat(sprintf("  Condition known      : %d / %d\n", sum(!is.na(master$Condition)), nrow(master)))
cat(sprintf("  Phenotype (curated)  : %d / %d\n", sum(!is.na(master$Phenotype)), nrow(master)))

cat("\n=== group sizes ===\n")
cat("Condition:\n");                  print(table(master$Condition, useNA = "ifany"))
cat("Phenotype (curated):\n");        print(table(master$Phenotype, useNA = "ifany"))
cat("Sex x Condition:\n");            print(table(master$Sex, master$Condition, useNA = "ifany"))

cat("\nwrote:", file.path(ENRICH, "sleap_all_batches_wide.tsv"), "\n")
