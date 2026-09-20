# ============================================================================
# 22_sleap_feedback_on_labels.R
#
# Does the SLEAPanalyzer work change the SUS/RES classification?
#
# It can only do so through ONE channel. The six classifier components are
# NOR D2, sucrose preference, delta CORT, body-weight development, adrenal
# weight and spleen weight. Five of those are not video measures at all, so
# nothing we did to SLEAPanalyzer can touch them. Only NOR D2 is derived from
# tracking -- and in the shipped classifier it came from the MANUAL BORIS
# scoring, not from SLEAP.
#
# So the question is precise: if the SLEAP-derived D2 (with the recalibrated
# 4 cm contact detector and the shape-matched object handling) is substituted
# for the BORIS D2, does any animal change label?
#
# Two things bound the answer in advance:
#   - SLEAP and BORIS agree closely on NOR (validated earlier: r 0.95-0.99,
#     CCC 0.93-0.97 on durations), so the substituted values are similar.
#   - NOR carries only 12.6% of the male composite and 1.9% of the female
#     one, so even a large change in NOR moves the composite little.
# Both are checked here rather than assumed.
# ============================================================================

suppressMessages({library(dplyr); library(tidyr); library(readxl)})

XL   <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis/E9_Behavior_Data.xlsx"
ANA  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis"
PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
RES  <- file.path(PROJ, "results")
SIX  <- c("NOR", "SucPref", "dCORT", "Bodyweight", "Adrenal", "Spleen")
canon <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(grepl("^[0-9]+$", x), sub("^0+(?=[0-9])", "", x, perl = TRUE), x)
}
sdp <- function(v) { v <- v[!is.na(v)]; sqrt(mean((v - mean(v))^2)) }

nb <- suppressMessages(read_excel(XL, sheet = "DLSsingleSlim_noBatch")) %>%
  mutate(across(-any_of(c("ID", "Group", "Sex", "Batch")), ~suppressWarnings(as.numeric(.))),
         key = canon(ID), Batch = paste0("B", Batch))
sus <- canon(trimws(readLines(file.path(ANA, "sus_animals.txt"), warn = FALSE)))
sl <- read.delim(file.path(PROJ, "enriched", "sleap_all_batches_wide.tsv"),
                 stringsAsFactors = FALSE, na.strings = "") %>%
  mutate(key = canon(ID)) %>% select(key, sleap_NOR_D2, sleap_NOR_qcPass)

d <- nb %>% left_join(sl, by = "key")

# ---- 1. coverage -----------------------------------------------------------
cat("======== 1. is there a SLEAP D2 for every animal? ========\n")
cat(sprintf("  animals in the classifier table : %d\n", nrow(d)))
cat(sprintf("  with a SLEAP NOR D2             : %d\n", sum(!is.na(d$sleap_NOR_D2))))
cat(sprintf("  missing                         : %d\n", sum(is.na(d$sleap_NOR_D2))))
if (any(is.na(d$sleap_NOR_D2)))
  print(d %>% filter(is.na(sleap_NOR_D2)) %>% select(ID, Batch, Sex, Group) %>%
          as.data.frame(), row.names = FALSE)

# ---- 2. how well do the two D2 measures agree? -----------------------------
cat("\n======== 2. BORIS D2 (in the classifier) vs SLEAP D2 ========\n")
cat("The stored NOR column is a z-score of the BORIS D2, so compare on ranks\n")
cat("and on the z-scale after applying the SAME transform to the SLEAP D2.\n\n")
ok <- d %>% filter(!is.na(sleap_NOR_D2), !is.na(NOR))
cat(sprintf("  Pearson r (stored NOR z vs raw SLEAP D2) : %.3f  (n = %d)\n",
            cor(ok$NOR, ok$sleap_NOR_D2), nrow(ok)))
cat(sprintf("  Spearman rho                             : %.3f\n",
            cor(ok$NOR, ok$sleap_NOR_D2, method = "spearman")))
cat("  positive r confirms the classifier orients D2 so higher = better.\n")

# ---- 3. build the SLEAP-substituted component ------------------------------
# Same rule the workbook used: z against that SEX's controls, population SD.
d$NOR_sleap_z <- NA_real_
for (sx in unique(d$Sex)) {
  i <- d$Sex == sx
  cv <- d$sleap_NOR_D2[i & d$Group == "CON"]
  d$NOR_sleap_z[i] <- (d$sleap_NOR_D2[i] - mean(cv, na.rm = TRUE)) / sdp(cv)
}
cat(sprintf("\n  agreement of the two z-scored NOR components: r = %.3f\n",
            cor(d$NOR, d$NOR_sleap_z, use = "complete.obs")))

# ---- 4. relabel with SLEAP NOR ---------------------------------------------
sis <- d %>% filter(Group == "SIS")
sis$A <- ifelse(sis$key %in% sus, "SUS", "RES")
NS <- tapply(sis$A == "SUS", sis$Sex, sum)
rel <- function(S, z) {
  o <- rep(NA_character_, nrow(S))
  for (sx in unique(S$Sex)) {
    i <- which(S$Sex == sx)
    o[i] <- ifelse(rank(z[i], ties.method = "first") <= NS[[sx]], "SUS", "RES")
  }
  o
}
sis$comp_boris <- rowMeans(sis[, SIX], na.rm = TRUE)
alt <- sis; alt$NOR <- ifelse(is.na(alt$NOR_sleap_z), alt$NOR, alt$NOR_sleap_z)
sis$comp_sleap <- rowMeans(alt[, SIX], na.rm = TRUE)
sis$L_sleap <- rel(sis, sis$comp_sleap)

cat("\n======== 3. does substituting the SLEAP D2 change any label? ========\n")
cat(sprintf("  composite correlation, BORIS-NOR vs SLEAP-NOR : r = %.4f\n",
            cor(sis$comp_boris, sis$comp_sleap, use = "complete.obs")))
cat(sprintf("  labels changed (SUS count held per sex)       : %d of %d\n",
            sum(sis$L_sleap != sis$A), nrow(sis)))
ch <- sis %>% filter(L_sleap != A) %>%
  mutate(across(c(comp_boris, comp_sleap, NOR, NOR_sleap_z), ~round(., 3))) %>%
  select(ID, Batch, Sex, current = A, with_sleap = L_sleap,
         NOR_boris = NOR, NOR_sleap = NOR_sleap_z, comp_boris, comp_sleap)
if (nrow(ch)) print(as.data.frame(ch), row.names = FALSE) else
  cat("  (none -- the classification is unchanged)\n")

# ---- 5. sensitivity: how much WOULD NOR have to move? ----------------------
cat("\n======== 4. how sensitive is the label to NOR at all? ========\n")
cat("Dropping NOR entirely is the largest possible perturbation to it.\n")
FIVE <- setdiff(SIX, "NOR")
sis$comp_noNOR <- rowMeans(sis[, FIVE], na.rm = TRUE)
sis$L_noNOR <- rel(sis, sis$comp_noNOR)
cat(sprintf("  labels changed if NOR is removed completely : %d of %d\n",
            sum(sis$L_noNOR != sis$A), nrow(sis)))
print(sis %>% group_by(Sex) %>%
        summarise(n = n(), changed_with_sleap = sum(L_sleap != A),
                  changed_if_NOR_dropped = sum(L_noNOR != A), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

write.csv(sis %>% select(ID, key, Batch, Sex, A, L_sleap, L_noNOR,
                         NOR, NOR_sleap_z, comp_boris, comp_sleap),
          file.path(RES, "sleap_feedback_on_labels.csv"), row.names = FALSE)
cat("\nwrote:", file.path(RES, "sleap_feedback_on_labels.csv"), "\n")
