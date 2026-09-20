# ============================================================================
# 19_verify_classification.R
#
# Re-derives the susceptibility classification from first principles and
# checks every assumption, because two earlier scripts got the base wrong.
#
# WHAT WAS WRONG BEFORE
#
# 18_within_sex_composite.R built labels on DLSsingleSlim (batch z-scored)
# and compared them with sus_animals.txt, which is NOT from that sheet. It
# then reported that "within-sex standardisation" moved 20 animals. Both the
# base and the premise were wrong:
#
#   - sus_animals.txt comes from DLSsingleSlim_noBatch, and that sheet is
#     ALREADY z-scored within sex (each sex has 12 controls, mean 0, sample
#     SD sqrt(12/11) = 1.0445).
#   - the list also uses a SEPARATE CUTOFF PER SEX. With those two facts the
#     list reproduces 93/93 exactly.
#
# So the classification was already within-sex on centring, on scaling and on
# the threshold. The "20 animals change" figure was an artefact of scoring on
# the wrong sheet.
#
# WHAT IS ACTUALLY STILL OPEN, and is tested here:
#   S1  the scale reference is the 12 CONTROLS per sex. Among the STRESSED
#       animals the components may still spread unequally, so a component can
#       still dominate the composite more in one sex than in the other.
#   S2  batch is not corrected. Each sex spans 3 batches and the controls of
#       those batches differ (F = 4.91, p = 0.0052), so within-sex batch
#       shifts survive into the score.
#   S3  the two cutoffs are not the same severity (f -0.227, m -0.444), so
#       "susceptible" is a more lenient label in females than in males.
# ============================================================================

suppressMessages({library(dplyr); library(tidyr); library(readxl)})

XL   <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis/E9_Behavior_Data.xlsx"
ANA  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis"
RES  <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris/results"
SIX  <- c("NOR", "SucPref", "dCORT", "Bodyweight", "Adrenal", "Spleen")
canon <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(grepl("^[0-9]+$", x), sub("^0+(?=[0-9])", "", x, perl = TRUE), x)
}

nb <- suppressMessages(read_excel(XL, sheet = "DLSsingleSlim_noBatch")) %>%
  mutate(across(-any_of(c("ID", "Group", "Sex", "Batch")), ~suppressWarnings(as.numeric(.))),
         key = canon(ID), Batch = paste0("B", Batch))
sus <- canon(trimws(readLines(file.path(ANA, "sus_animals.txt"), warn = FALSE)))

# ---- V. exact reproduction of the shipped list -----------------------------
sis <- nb %>% filter(Group == "SIS")
sis$comp  <- rowMeans(sis[, SIX], na.rm = TRUE)
sis$yours <- ifelse(sis$key %in% sus, "SUS", "RES")

cat("======== V. reproducing sus_animals.txt exactly ========\n")
cuts <- list()
for (sx in c("f", "m")) {
  i <- sis$Sex == sx
  k <- sum(sis$yours[i] == "SUS")
  v <- sort(sis$comp[i])
  cuts[[sx]] <- c(lo = v[k], hi = v[k + 1])
  cat(sprintf("  %s: %2d SUS of %2d. any cutoff in (%.4f, %.4f) gives your exact labels\n",
              sx, k, sum(i), v[k], v[k + 1]))
}
sis$recon <- vapply(seq_len(nrow(sis)), function(r)
  if (sis$comp[r] <= mean(cuts[[sis$Sex[r]]])) "SUS" else "RES", character(1))
cat(sprintf("  reproduction: %d of %d (%.1f%%)\n",
            sum(sis$recon == sis$yours), nrow(sis), 100 * mean(sis$recon == sis$yours)))
cat("  => your list = mean of the six sex-z-scored components, thresholded within sex.\n")

# ---- S1. what each component actually contributes --------------------------
# C = (1/k) * sum(z_i), so var(C) = (1/k) * sum cov(z_i, C) and a component's
# share of the composite is cov(z_i, C) / (k * var(C)). These sum to 1.
# Two ways to get this wrong, both of which I did earlier:
#   - var(z_i)/sum(var(z_j)) ignores the covariance terms entirely, and is
#     tautologically 1/k once the components are scaled to unit variance.
#   - cov(z_i, C)/var(C) without the factor k sums to k, not 1.
# The assertion below makes either failure loud instead of silent.
influence <- function(df, cols) {
  M <- as.matrix(df[, cols])
  M <- M[stats::complete.cases(M), , drop = FALSE]
  C <- rowMeans(M)
  sh <- 100 * sapply(cols, function(cc) cov(M[, cc], C)) / (length(cols) * var(C))
  stopifnot(abs(sum(sh) - 100) < 1e-6)
  round(sh, 1)
}
cat("\n======== S1. share of the composite each component drives (%) ========\n")
print(data.frame(component = SIX,
                 male   = influence(sis %>% filter(Sex == "m"), SIX),
                 female = influence(sis %>% filter(Sex == "f"), SIX)), row.names = FALSE)
cat("equal weight would be 16.7%. Shares sum to 100 by construction.\n")

# ---- alternative scores ----------------------------------------------------
# E1: additionally divide each component by its within-sex SD among the
#     STRESSED animals, so each contributes equally in each sex.
e1 <- sis
sdv <- sis %>% group_by(Sex) %>%
  summarise(across(all_of(SIX), ~sd(., na.rm = TRUE)), .groups = "drop")
for (sx in unique(sis$Sex)) {
  i <- e1$Sex == sx
  for (cc in SIX) e1[[cc]][i] <- e1[[cc]][i] / sdv[[cc]][sdv$Sex == sx]
}

# E2: additionally centre each component on its OWN BATCH's controls,
#     removing within-sex batch shifts. Units stay the sex-control SD.
bmean <- nb %>% filter(Group == "CON") %>% group_by(Batch) %>%
  summarise(across(all_of(SIX), ~mean(., na.rm = TRUE)), .groups = "drop")
e2 <- sis
for (b in bmean$Batch) {
  i <- e2$Batch == b
  for (cc in SIX) e2[[cc]][i] <- e2[[cc]][i] - bmean[[cc]][bmean$Batch == b]
}

# E3: both.
e3 <- e2
sdv3 <- e2 %>% group_by(Sex) %>%
  summarise(across(all_of(SIX), ~sd(., na.rm = TRUE)), .groups = "drop")
for (sx in unique(e3$Sex)) {
  i <- e3$Sex == sx
  for (cc in SIX) e3[[cc]][i] <- e3[[cc]][i] / sdv3[[cc]][sdv3$Sex == sx]
}

cat("\nafter equalising component spread within sex (E1):\n")
print(data.frame(component = SIX,
                 male   = influence(e1 %>% filter(Sex == "m"), SIX),
                 female = influence(e1 %>% filter(Sex == "f"), SIX)), row.names = FALSE)

# ---- relabel, count-matched WITHIN SEX -------------------------------------
# Each sex keeps its own SUS count (f 22, m 17), so the comparison isolates
# the change in SCORE and never confounds it with a threshold choice.
relabel <- function(df) {
  z <- rowMeans(df[, SIX], na.rm = TRUE)
  out <- rep(NA_character_, nrow(df))
  for (sx in unique(df$Sex)) {
    i <- which(df$Sex == sx)
    k <- sum(sis$yours[i] == "SUS")
    out[i] <- ifelse(rank(z[i], ties.method = "first") <= k, "SUS", "RES")
  }
  out
}
sis$E1 <- relabel(e1); sis$E2 <- relabel(e2); sis$E3 <- relabel(e3)
sis$z_yours <- sis$comp
sis$z_E3 <- rowMeans(e3[, SIX], na.rm = TRUE)

cat("\n======== how many animals move, count-matched within sex ========\n")
lbl <- c(E1 = "E1  + equalise component spread within sex",
         E2 = "E2  + centre on each batch's own controls",
         E3 = "E3  both")
for (nm in c("E1", "E2", "E3"))
  cat(sprintf("%-46s %2d of 93 change (f %d, m %d)\n", lbl[nm],
              sum(sis[[nm]] != sis$yours),
              sum(sis[[nm]] != sis$yours & sis$Sex == "f"),
              sum(sis[[nm]] != sis$yours & sis$Sex == "m")))

cat("\n======== animals classified differently from your list ========\n")
ch <- sis %>%
  filter(E3 != yours | E1 != yours | E2 != yours) %>%
  mutate(across(c(z_yours, z_E3), ~round(., 3))) %>%
  select(ID, Batch, Sex, yours, E1, E2, E3, z_yours, z_E3) %>%
  arrange(Sex, yours, ID)
print(as.data.frame(ch), row.names = FALSE)

# ---- S3. the severity asymmetry --------------------------------------------
cat("\n======== S3. is 'susceptible' the same severity in both sexes? ========\n")
cat(sprintf("  female cutoff ~ %.3f -> %.1f%% of stressed females called SUS\n",
            mean(cuts$f), 100 * mean(sis$yours[sis$Sex == "f"] == "SUS")))
cat(sprintf("  male   cutoff ~ %.3f -> %.1f%% of stressed males   called SUS\n",
            mean(cuts$m), 100 * mean(sis$yours[sis$Sex == "m"] == "SUS")))
print(sis %>% group_by(Sex) %>%
        summarise(n = n(), SUS = sum(yours == "SUS"),
                  pct = round(100 * mean(yours == "SUS"), 1),
                  median_comp = round(median(comp), 3), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

write.csv(sis %>% select(ID, key, Batch, Sex, z_yours, z_E3, yours, E1, E2, E3),
          file.path(RES, "classification_verified.csv"), row.names = FALSE)
cat("\nwrote:", file.path(RES, "classification_verified.csv"), "\n")
