# ============================================================================
# 16_reproduce_classifier.R
#
# The actual scheme, per the experimenter: each batch was CENTRED on its own
# controls, but the SD used to scale was POOLED ACROSS ALL BATCHES -- not the
# 4-animal per-batch SD I assumed in 15_phenotype_schemes.R. And the previous
# outcome analyses did NOT include batch as a covariate.
#
# "SD pooled across all batches" is ambiguous in one important way, so both
# readings are tested:
#
#   pooled_raw       sd() of all control values on the raw scale. Batch means
#                    differ, so this absorbs between-batch spread INTO the
#                    scale and shrinks every z-score.
#   pooled_within    controls centred within batch first, then pooled. This is
#                    the textbook pooled SD: between-batch shifts are removed
#                    before estimating spread, and it has 24 - 6 = 18 df.
#
# Whichever reproduces sus_animals_batchCorrected is the one that was used.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
})

PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
SIS  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis"
ANA  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis"
RES  <- file.path(PROJ, "results")

canon <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(grepl("^[0-9]+$", x), sub("^0+(?=[0-9])", "", x, perl = TRUE), x)
}
SIX  <- c("NOR", "SucPref", "dCORT", "Bodyweight", "Adrenal", "Spleen")
FIVE <- setdiff(SIX, "NOR")

comp <- suppressMessages(read_excel(file.path(SIS, "E9_Behavior_Data.xlsx"),
                                    sheet = "DLSsingleSlim")) %>%
  mutate(across(-c(ID, Group, Sex, Batch), ~suppressWarnings(as.numeric(.))),
         key = canon(ID), Batch = paste0("B", Batch))

sus_orig <- canon(trimws(readLines(file.path(ANA, "sus_animals.txt"), warn = FALSE)))
sus_bc   <- canon(trimws(readLines(file.path(ANA, "sus_animals_batchCorrected.txt"), warn = FALSE)))

# --- scaling variants -------------------------------------------------------
# centre: "batch" = each batch's own control mean; "none" = no centring at all
#         (the uncorrected list)
# scale : per-batch SD | pooled raw | pooled within-batch | within-sex pooled
zscore <- function(d, col, centre = "batch", scale = "pooled_within") {
  x <- d[[col]]
  is_con <- d$Group == "CON"
  out <- rep(NA_real_, nrow(d))

  s_global <- switch(scale,
    pooled_raw = sd(x[is_con & is.finite(x)]),
    pooled_within = {
      r <- unlist(lapply(split(which(is_con), d$Batch[is_con]), function(i) {
        v <- x[i]; v <- v[is.finite(v)]; if (length(v) < 2) NULL else v - mean(v)
      }))
      # pooled within-batch SD: sum of squared deviations / (N - n_batches)
      nb <- length(unique(d$Batch[is_con]))
      sqrt(sum(r^2) / (length(r) - nb))
    },
    NA_real_)

  for (b in unique(d$Batch)) {
    ib <- d$Batch == b
    cv <- x[ib & is_con]; cv <- cv[is.finite(cv)]
    if (length(cv) < 2) next
    m <- if (centre == "batch") mean(cv) else mean(x[is_con & is.finite(x)])
    s <- switch(scale,
      per_batch = sd(cv),
      within_sex = { sx <- d$Sex[ib][1]
                     v <- x[is_con & d$Sex == sx]; sd(v[is.finite(v)]) },
      s_global)
    if (!is.finite(s) || s == 0) next
    out[ib] <- (x[ib] - m) / s
  }
  out
}

build <- function(cols, centre, scale)
  rowMeans(vapply(cols, function(cc) zscore(comp, cc, centre, scale),
                  numeric(nrow(comp))), na.rm = TRUE)

VARIANTS <- list(
  "no centring, pooled raw SD"        = list(SIX, "none",  "pooled_raw"),
  "batch-centred, pooled RAW SD"      = list(SIX, "batch", "pooled_raw"),
  "batch-centred, pooled WITHIN SD"   = list(SIX, "batch", "pooled_within"),
  "batch-centred, within-sex SD"      = list(SIX, "batch", "within_sex"),
  "batch-centred, per-batch SD"       = list(SIX, "batch", "per_batch")
)
for (nm in names(VARIANTS)) {
  v <- VARIANTS[[nm]]
  comp[[nm]] <- build(v[[1]], v[[2]], v[[3]])
}

sis <- comp %>% filter(Group == "SIS")
sis$L_orig <- ifelse(sis$key %in% sus_orig, "SUS", "RES")
sis$L_bc   <- ifelse(sis$key %in% sus_bc,   "SUS", "RES")

# --- which variant reproduces the existing lists? ---------------------------
auc <- function(score, lab) {
  ok <- is.finite(score); score <- score[ok]; lab <- lab[ok]
  r <- rank(score); n1 <- sum(lab); n0 <- sum(!lab)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  1 - (sum(r[lab]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
best_acc <- function(score, lab) {
  ok <- is.finite(score); score <- score[ok]; lab <- lab[ok]
  cand <- sort(unique(score))
  a <- vapply(cand, function(k) mean((score <= k) == lab), numeric(1))
  c(cut = cand[which.max(a)], acc = max(a), n_mis = round((1 - max(a)) * length(lab)))
}

cat("=== which scaling reproduces each shipped list? ===\n")
rep_tab <- bind_rows(lapply(names(VARIANTS), function(nm) {
  bo <- best_acc(sis[[nm]], sis$L_orig == "SUS")
  bb <- best_acc(sis[[nm]], sis$L_bc   == "SUS")
  data.frame(variant = nm,
             AUC_roster = auc(sis[[nm]], sis$L_orig == "SUS"),
             acc_roster = bo[["acc"]], mis_roster = bo[["n_mis"]],
             AUC_bc = auc(sis[[nm]], sis$L_bc == "SUS"),
             acc_bc = bb[["acc"]], mis_bc = bb[["n_mis"]], cut_bc = bb[["cut"]])
}))
print(rep_tab %>% mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(),
      row.names = FALSE)
cat("\nmis_* = animals the variant cannot place on the right side of any cutoff.\n")

# --- the reproducing variant, then leave-NOR-out ---------------------------
best_v <- rep_tab$variant[which.max(rep_tab$acc_bc)]
cat(sprintf("\nbest reproduction of sus_animals_batchCorrected: %s (%.1f%%, %d misplaced)\n",
            best_v, 100 * max(rep_tab$acc_bc), rep_tab$mis_bc[which.max(rep_tab$acc_bc)]))

spec <- VARIANTS[[best_v]]
CUT <- rep_tab$cut_bc[which.max(rep_tab$acc_bc)]
sis$z6  <- build(SIX,  spec[[2]], spec[[3]])[comp$Group == "SIS"]
sis$z5  <- build(FIVE, spec[[2]], spec[[3]])[comp$Group == "SIS"]
sis$L_reproduced <- ifelse(sis$z6 <= CUT, "SUS", "RES")
sis$L_loo        <- ifelse(sis$z5 <= CUT, "SUS", "RES")

cat(sprintf("\nreproduced label: %d SUS (shipped bc list: %d) | leave-NOR-out: %d SUS\n",
            sum(sis$L_reproduced == "SUS"), sum(sis$L_bc == "SUS"),
            sum(sis$L_loo == "SUS")))
cat(sprintf("leave-NOR-out differs from the reproduced label for %d of %d animals\n",
            sum(sis$L_loo != sis$L_reproduced), nrow(sis)))

# --- outcome: the 2x2 the experimenter asked for ---------------------------
dat <- read.delim(file.path(PROJ, "enriched", "sleap_all_batches_wide.tsv"),
                  stringsAsFactors = FALSE, na.strings = "") %>% mutate(key = canon(ID))
j <- dat %>% inner_join(sis %>% select(key, L_orig, L_bc, L_reproduced, L_loo), by = "key")

fit_one <- function(lc, metric, use_batch, tag) {
  d <- j %>% filter(!is.na(.data[[metric]]), .data[[lc]] %in% c("RES", "SUS"))
  if (nrow(d) < 10) return(NULL)
  d$g <- factor(d[[lc]], levels = c("RES", "SUS"))
  fit <- lm(reformulate(if (use_batch) c("g", "Batch") else "g", metric), data = d)
  co <- summary(fit)$coefficients["gSUS", ]; ci <- confint(fit)["gSUS", ]
  s <- summary(fit)$sigma
  data.frame(label = tag, batch_in_model = use_batch, metric = metric, n = nrow(d),
             n_sus = sum(d$g == "SUS"), d = co[1] / s,
             lo = ci[1] / s, hi = ci[2] / s, p = co[4], stringsAsFactors = FALSE)
}

LAB <- c(L_orig = "uncorrected roster (NOR inside)",
         L_bc = "batchCorrected, shipped (NOR inside)",
         L_loo = "reproduced scheme, leave-NOR-out")
out <- bind_rows(lapply(names(LAB), function(lc)
  bind_rows(lapply(c(TRUE, FALSE), function(ub)
    fit_one(lc, "sleap_NOR_D2", ub, LAB[[lc]])))))
write.csv(out, file.path(RES, "classifier_reproduction_outcomes.csv"), row.names = FALSE)

cat("\n=== NOR D2 under each label, with and without batch in the model ===\n")
print(out %>% mutate(d = round(d, 2), CI = sprintf("[%.2f, %.2f]", lo, hi),
                     p = signif(p, 3)) %>%
        select(label, batch_in_model, n_sus, d, CI, p) %>% as.data.frame(),
      row.names = FALSE)
cat("\nOnly the leave-NOR-out rows are interpretable; the others contain NOR.\n")

write.csv(sis %>% select(ID, key, Batch, Sex, z6, z5, L_orig, L_bc, L_reproduced, L_loo),
          file.path(RES, "classifier_reproduction_labels.csv"), row.names = FALSE)
cat("\nwrote:", file.path(RES, "classifier_reproduction_labels.csv"), "\n")
