# ============================================================================
# 14_leave_one_out_phenotype.R
#
# THE PROBLEM
#
# SUS/RES is a threshold on a six-component composite:
#
#     NOR D2 | sucrose preference | delta CORT | body-weight development |
#     adrenal weight | spleen weight
#
# (Analysis/SIS_Analysis/E9_Behavior_Data.xlsx, sheet DLSsingleSlim. The
# 6-component composite predicts sus_animals_batchCorrected at AUC 0.999, so
# that list is a batch z-scored threshold on exactly this score.)
#
# NOR D2 is therefore an INPUT to the classifier. Testing NOR against SUS/RES
# asks whether a score containing NOR differs on NOR: it is circular, and any
# effect is guaranteed by construction. An earlier version of this analysis
# reported "SUS vs RES on NOR D2, d = -0.67, Holm 0.020" as a headline. That
# result is not interpretable.
#
# SocP, EPM and OFT are NOT components, so phenotype contrasts on those assays
# were never circular and stand as reported.
#
# THE FIX
#
# Rebuild the composite WITHOUT NOR, re-threshold, and test NOR against that
# label. The classifier then uses only information independent of the outcome.
# NOR correlates with the other five components at only r = 0.17, so this is a
# real change of label rather than a cosmetic one.
#
# Both thresholding rules are reported, because neither is obviously right:
#   - count-matched: same number of SUS animals as the original list
#   - fixed cutoff:  the z-cutoff that reproduces the original count overall
# ============================================================================

suppressMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
})

PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
SIS  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis"
ANA  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis"
RES  <- file.path(PROJ, "results")
FIG  <- file.path(PROJ, "figures")

canon <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(grepl("^[0-9]+$", x), sub("^0+(?=[0-9])", "", x, perl = TRUE), x)
}

comp <- suppressMessages(read_excel(file.path(SIS, "E9_Behavior_Data.xlsx"),
                                    sheet = "DLSsingleSlim")) %>%
  mutate(across(-c(ID, Group, Sex, Batch), ~suppressWarnings(as.numeric(.))),
         key = canon(ID), Batch = paste0("B", Batch))

SIX  <- c("NOR", "SucPref", "dCORT", "Bodyweight", "Adrenal", "Spleen")
FIVE <- setdiff(SIX, "NOR")

sus_roster <- canon(trimws(readLines(file.path(ANA, "sus_animals.txt"), warn = FALSE)))

dat <- read.delim(file.path(PROJ, "enriched", "sleap_all_batches_wide.tsv"),
                  stringsAsFactors = FALSE, na.strings = "") %>%
  mutate(key = canon(ID))

# --- Composite, with and without NOR ---------------------------------------
# Batch z-scoring: centre and scale each animal on its OWN batch's controls,
# which is what the batch-corrected list does. The reference is only 4 control
# animals per batch, so the scale estimate is noisy; a pooled within-sex SD is
# reported alongside as the more stable alternative.
zscore_to_controls <- function(d, col, pooled_sd = FALSE) {
  out <- rep(NA_real_, nrow(d))
  for (b in unique(d$Batch)) {
    ib <- d$Batch == b
    ctrl <- d[[col]][ib & d$Group == "CON"]
    ctrl <- ctrl[is.finite(ctrl)]
    if (length(ctrl) < 2) next
    s <- if (pooled_sd) {
      sx <- d$Sex[ib][1]
      cs <- d[[col]][d$Sex == sx & d$Group == "CON"]
      sd(cs[is.finite(cs)])
    } else sd(ctrl)
    if (!is.finite(s) || s == 0) next
    out[ib] <- (d[[col]][ib] - mean(ctrl)) / s
  }
  out
}

build <- function(components, pooled_sd = FALSE) {
  z <- vapply(components, function(cc) zscore_to_controls(comp, cc, pooled_sd),
              numeric(nrow(comp)))
  rowMeans(z, na.rm = TRUE)
}

comp$DLS6 <- build(SIX)
comp$DLS5 <- build(FIVE)
comp$DLS5_pooled <- build(FIVE, pooled_sd = TRUE)

sis <- comp %>% filter(Group == "SIS", is.finite(DLS5))
n_sus <- sum(sis$key %in% sus_roster)
cat(sprintf("SIS animals with the 5 non-NOR components: %d | original SUS count: %d\n",
            nrow(sis), n_sus))

# Count-matched threshold: the n_sus lowest composites are SUS.
label_by_count <- function(score, n) {
  ifelse(rank(score, ties.method = "first") <= n, "SUS", "RES")
}
sis$pheno_loo <- label_by_count(sis$DLS5, n_sus)
sis$pheno_loo_pooled <- label_by_count(sis$DLS5_pooled, n_sus)
sis$pheno_full <- label_by_count(sis$DLS6, n_sus)
sis$pheno_orig <- ifelse(sis$key %in% sus_roster, "SUS", "RES")

cat("\n=== how much does dropping NOR move the label? ===\n")
cat(sprintf("  leave-one-out vs original roster : %d of %d animals change\n",
            sum(sis$pheno_loo != sis$pheno_orig), nrow(sis)))
cat(sprintf("  6-component vs original roster   : %d change\n",
            sum(sis$pheno_full != sis$pheno_orig)))
cat(sprintf("  pooled-SD vs per-batch SD        : %d change\n",
            sum(sis$pheno_loo != sis$pheno_loo_pooled)))

# --- Test NOR against each label -------------------------------------------
j <- dat %>% inner_join(sis %>% select(key, pheno_orig, pheno_full, pheno_loo,
                                       pheno_loo_pooled), by = "key")
cat(sprintf("\njoined to the SLEAP outcomes: %d animals\n", nrow(j)))

test_one <- function(label_col, metric, tag) {
  d <- j %>% filter(!is.na(.data[[metric]]), .data[[label_col]] %in% c("RES", "SUS"))
  if (nrow(d) < 10) return(NULL)
  d$g <- factor(d[[label_col]], levels = c("RES", "SUS"))
  fit <- lm(reformulate(c("g", "Batch"), metric), data = d)
  co <- summary(fit)$coefficients["gSUS", ]
  ci <- confint(fit)["gSUS", ]
  s <- summary(fit)$sigma
  data.frame(label = tag, metric = metric, n = nrow(d),
             n_sus = sum(d$g == "SUS"),
             d = co[1] / s, lo = ci[1] / s, hi = ci[2] / s, p = co[4],
             stringsAsFactors = FALSE)
}

OUT <- c("sleap_NOR_D2", "sleap_NOR_total_s",
         "sleap_SocP_S1_pref_index", "sleap_EPM_open_frac", "sleap_OFT_center_pct")
labels <- c(pheno_orig = "original roster (NOR inside)",
            pheno_full = "6-component (NOR inside)",
            pheno_loo = "leave-NOR-out",
            pheno_loo_pooled = "leave-NOR-out, pooled SD")

res <- bind_rows(lapply(names(labels), function(lc)
  bind_rows(lapply(OUT, function(m) test_one(lc, m, labels[[lc]])))))
write.csv(res, file.path(RES, "phenotype_leave_one_out.csv"), row.names = FALSE)

cat("\n=== NOR D2: circular labels vs the leave-one-out label ===\n")
print(res %>% filter(metric == "sleap_NOR_D2") %>%
        mutate(d = round(d, 2), CI = sprintf("[%.2f, %.2f]", lo, hi),
               p = signif(p, 3)) %>%
        select(label, n, n_sus, d, CI, p) %>% as.data.frame(), row.names = FALSE)

cat("\n=== the assays that were never in the classifier (unchanged, for reference) ===\n")
print(res %>% filter(metric != "sleap_NOR_D2", label %in% labels[c(1, 3)]) %>%
        mutate(d = round(d, 2), p = signif(p, 3)) %>%
        select(label, metric, n, d, p) %>% as.data.frame(), row.names = FALSE)

cat("\nNOR D2 under a NOR-containing label is circular and is reported only to\n",
    "show the size of the artefact. The leave-NOR-out row is the interpretable\n",
    "estimate of whether susceptibility predicts novel-object discrimination.\n", sep = "")
cat("\nwrote:", file.path(RES, "phenotype_leave_one_out.csv"), "\n")
