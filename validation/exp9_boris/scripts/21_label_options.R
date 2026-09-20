# ============================================================================
# 21_label_options.R
#
# The decision: keep sus_animals.txt, or adopt one of the corrected variants?
#
# THE OPTIONS  (all six components, all thresholded within sex)
#
#   A  YOURS   centre + scale on that SEX's 12 controls. Reproduces
#              sus_animals.txt 93/93. Already within-sex on every axis.
#   E1 A + divide each component by its within-sex SD among the STRESSED
#              animals, so each contributes equally to the composite.
#   E2 A + centre on that BATCH's own 4 controls instead of the sex pool,
#              removing within-sex batch shifts in the reference.
#   E3 both.
#
# Effect size is not a valid criterion for choosing between these -- picking
# the labelling that produces the nicest p-values is circular. Two criteria
# that ARE valid:
#
#   BIAS      does the scheme remove a known artefact? E2 does: under A a
#             batch whose controls happened to perform well pushes its own
#             stressed animals toward SUS.
#   VARIANCE  how stable is the label? Every scheme estimates its reference
#             from a handful of controls -- 12 per sex under A, only 4 per
#             batch under E2. A scheme that fixes a bias by spending a lot of
#             variance is not obviously an improvement.
#
# The bootstrap below measures VARIANCE directly: resample the control
# animals that define each reference, rebuild the score, relabel, and count
# how often each animal flips. That is the number that should decide this.
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
con <- nb %>% filter(Group == "CON")
sis <- nb %>% filter(Group == "SIS")
sis$yours <- ifelse(sis$key %in% sus, "SUS", "RES")
NSUS <- tapply(sis$yours == "SUS", sis$Sex, sum)   # f 22, m 17

# --- score builders ---------------------------------------------------------
# Each returns the composite for the SIS animals, given a control table.
# All work on the _noBatch scale; since that is an affine transform of the raw
# values, recentring/rescaling here is equivalent to doing it on raw data.
score <- function(S, C, batch_centre, equalise) {
  M <- as.matrix(S[, SIX])
  for (cc in seq_along(SIX)) {
    v <- SIX[cc]
    for (sx in unique(S$Sex)) {
      i <- S$Sex == sx
      s_ref <- sd(C[[v]][C$Sex == sx], na.rm = TRUE)
      if (!is.finite(s_ref) || s_ref == 0) s_ref <- 1
      if (batch_centre) {
        for (b in unique(S$Batch[i])) {
          j <- i & S$Batch == b
          m <- mean(C[[v]][C$Batch == b], na.rm = TRUE)
          if (!is.finite(m)) m <- mean(C[[v]][C$Sex == sx], na.rm = TRUE)
          M[j, cc] <- (M[j, cc] - m) / s_ref
        }
      } else {
        m <- mean(C[[v]][C$Sex == sx], na.rm = TRUE)
        M[i, cc] <- (M[i, cc] - m) / s_ref
      }
    }
  }
  if (equalise) for (sx in unique(S$Sex)) {
    i <- S$Sex == sx
    for (cc in seq_along(SIX)) {
      s2 <- sd(M[i, cc], na.rm = TRUE)
      if (is.finite(s2) && s2 > 0) M[i, cc] <- M[i, cc] / s2
    }
  }
  rowMeans(M, na.rm = TRUE)
}

relabel <- function(S, z) {
  out <- rep(NA_character_, nrow(S))
  for (sx in unique(S$Sex)) {
    i <- which(S$Sex == sx)
    out[i] <- ifelse(rank(z[i], ties.method = "first") <= NSUS[[sx]], "SUS", "RES")
  }
  out
}

OPT <- list(A  = c(batch = FALSE, eq = FALSE),
            E1 = c(batch = FALSE, eq = TRUE),
            E2 = c(batch = TRUE,  eq = FALSE),
            E3 = c(batch = TRUE,  eq = TRUE))
Z <- L <- list()
for (o in names(OPT)) {
  Z[[o]] <- score(sis, con, OPT[[o]][["batch"]], OPT[[o]][["eq"]])
  L[[o]] <- relabel(sis, Z[[o]])
}
stopifnot(identical(L$A, sis$yours))   # option A must BE your list
cat("check: option A reproduces sus_animals.txt exactly -> OK\n")

# --- 1. agreement matrix ----------------------------------------------------
cat("\n======== 1. how much do the options disagree? ========\n")
ag <- outer(names(OPT), names(OPT),
            Vectorize(function(a, b) sum(L[[a]] != L[[b]])))
dimnames(ag) <- list(names(OPT), names(OPT))
cat("animals labelled differently (of 93):\n"); print(ag)

# --- 2. margin: how close is each animal to its own cutoff? -----------------
margin <- function(o) {
  z <- Z[[o]]; m <- rep(NA_real_, length(z))
  for (sx in unique(sis$Sex)) {
    i <- which(sis$Sex == sx); v <- sort(z[i]); k <- NSUS[[sx]]
    cut <- mean(c(v[k], v[k + 1]))
    m[i] <- (z[i] - cut) / sd(z[i])      # in within-sex SD of the composite
  }
  m
}
MA <- sapply(names(OPT), margin)
cat("\n======== 2. how many animals sit within noise of the cutoff? ========\n")
print(data.frame(option = names(OPT),
                 within_0.25_SD = colSums(abs(MA) < 0.25),
                 within_0.50_SD = colSums(abs(MA) < 0.50)), row.names = FALSE)
cat("The threshold is a cut through a continuum; animals this close to it are\n")
cat("not meaningfully 'susceptible' or 'resilient' under ANY scheme.\n")

# --- 3. bootstrap the control reference -------------------------------------
cat("\n======== 3. label stability: resampling the control reference ========\n")
set.seed(20260920)
NB <- 2000
flip <- matrix(0, nrow(sis), length(OPT), dimnames = list(NULL, names(OPT)))
for (b in seq_len(NB)) {
  # resample controls within the stratum that defines each reference
  idx_sex <- unlist(lapply(split(seq_len(nrow(con)), con$Sex),
                           function(i) sample(i, length(i), TRUE)))
  idx_bat <- unlist(lapply(split(seq_len(nrow(con)), con$Batch),
                           function(i) sample(i, length(i), TRUE)))
  for (o in names(OPT)) {
    C <- con[if (OPT[[o]][["batch"]]) idx_bat else idx_sex, ]
    lb <- relabel(sis, score(sis, C, OPT[[o]][["batch"]], OPT[[o]][["eq"]]))
    flip[, o] <- flip[, o] + (lb != L[[o]])
  }
}
flip <- flip / NB
cat(sprintf("%d bootstrap replicates.\n\n", NB))
print(data.frame(option = names(OPT),
                 mean_flip_rate = round(colMeans(flip), 3),
                 animals_over_20pct = colSums(flip > 0.20),
                 animals_over_33pct = colSums(flip > 0.33)), row.names = FALSE)
cat("\nmean_flip_rate = expected share of the 93 labels that change if the\n")
cat("control animals had been a different random draw from the same population.\n")

# --- 4. the animals that change --------------------------------------------
cat("\n======== 4. animals that change, with how solid each call is ========\n")
tab <- data.frame(ID = sis$ID, Batch = sis$Batch, Sex = sis$Sex,
                  A = L$A, E1 = L$E1, E2 = L$E2, E3 = L$E3,
                  margin_A = round(MA[, "A"], 2),
                  flip_A = round(flip[, "A"], 2),
                  flip_E2 = round(flip[, "E2"], 2),
                  stringsAsFactors = FALSE)
ch <- tab %>% filter(E1 != A | E2 != A | E3 != A) %>% arrange(Sex, A, ID)
print(as.data.frame(ch), row.names = FALSE)
cat("\nmargin_A = distance from your cutoff, in within-sex SD of the composite.\n")
cat("flip_*   = how often that animal changes label under resampling.\n")
cat(sprintf("\nof the %d animals that move under any option, %d already have\n",
            nrow(ch), sum(ch$flip_A > 0.20)))
cat("a >20%% flip rate under YOUR OWN scheme -- i.e. they were never firmly placed.\n")

write.csv(tab, file.path(RES, "label_options.csv"), row.names = FALSE)
cat("\nwrote:", file.path(RES, "label_options.csv"), "\n")
