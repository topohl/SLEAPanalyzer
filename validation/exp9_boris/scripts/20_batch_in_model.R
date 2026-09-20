# ============================================================================
# 20_batch_in_model.R
#
# Should Batch be a covariate in the OUTCOME model?
#
# This is a different question from whether to batch-correct the CLASSIFIER
# inputs (19_verify_classification.R, scheme E2). The two are independent
# choices and they are not double-correction: E2 changes WHO is called SUS,
# batch-in-the-model changes how SUS and RES are COMPARED on assays that
# played no part in the labelling.
#
# Adding Batch to the outcome model does two things at once:
#
#   (a) PRECISION. It removes between-batch variance from the residual. If
#       batch genuinely shifts the outcome, the SE shrinks and power rises.
#       This is always legitimate.
#   (b) ADJUSTMENT. It removes whatever part of the group effect is aligned
#       with batch. Whether that is desirable depends entirely on WHY label
#       and batch are aligned.
#
# The (b) question is the whole issue, and the answer differs by contrast:
#
#   CON vs SIS   Condition is the experimental assignment and is balanced
#                across batches by design. Batch cannot confound it, so
#                adding Batch is pure (a) -- free precision, no cost.
#
#   SUS vs RES   The label is DERIVED from data, so it can be entangled with
#                batch. Under the shipped list the reference is one sex-wide
#                control pool, so a batch whose controls happened to perform
#                well pushes its own stressed animals toward SUS. That
#                alignment is an artefact of the reference, not an excess of
#                real susceptibility -- so adjusting for it is appropriate.
#                Under E2 the labels are already batch-centred, so little
#                alignment is left and Batch reverts to being mostly (a).
#
# Sex is perfectly nested in batch (B1,B2,B5 male; B3,B4,B6 female), so in a
# pooled model Batch absorbs Sex entirely. A pooled model therefore cannot
# report a sex main effect; that is a property of the design, not a bug.
# ============================================================================

suppressMessages({library(dplyr); library(tidyr)})

PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
RES  <- file.path(PROJ, "results")
canon <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(grepl("^[0-9]+$", x), sub("^0+(?=[0-9])", "", x, perl = TRUE), x)
}

dat <- read.delim(file.path(PROJ, "enriched", "sleap_all_batches_wide.tsv"),
                  stringsAsFactors = FALSE, na.strings = "") %>%
  mutate(key = canon(ID))
lab <- read.csv(file.path(RES, "classification_verified.csv"),
                stringsAsFactors = FALSE) %>%
  mutate(key = canon(key)) %>% select(key, yours, E1, E2, E3)

d <- dat %>% left_join(lab, by = "key") %>%
  mutate(across(c(yours, E1, E2, E3), ~ifelse(Condition == "CON", "CON", .)),
         # the wide table says Female/Male, the label file says f/m
         Sex = recode(Sex, Female = "f", Male = "m"),
         cond = Condition)

# Outcomes that played NO part in the classification.
IND <- c("EPM open-arm fraction" = "sleap_EPM_open_frac",
         "EPM centre time"       = "sleap_EPM_center_s",
         "EPM nose dips"         = "sleap_EPM_nosedip_n",
         "EPM distance"          = "sleap_EPM_distance_cm",
         "OFT centre %"          = "sleap_OFT_center_pct",
         "OFT corner %"          = "sleap_OFT_corner_pct",
         "OFT distance"          = "sleap_OFT_distance_cm",
         "SocP S1 preference"    = "sleap_SocP_S1_pref_index",
         "SocP S2 preference"    = "sleap_SocP_S2_pref_index")

# ---------------------------------------------------------------------------
cat("======== 1. is the LABEL entangled with batch? ========\n")
cat("(within sex, because batch is nested in sex -- a pooled table would just\n")
cat(" be showing the sex difference again)\n\n")
cramer <- function(x, y) {
  tb <- table(x, y); if (any(dim(tb) < 2)) return(NA_real_)
  suppressWarnings(sqrt(chisq.test(tb)$statistic /
                          (sum(tb) * (min(dim(tb)) - 1)))[[1]])
}
ent <- bind_rows(lapply(c("yours", "E1", "E2", "E3"), function(s) {
  sis <- d %>% filter(Condition == "SIS", !is.na(.data[[s]]))
  bind_rows(lapply(c("f", "m"), function(sx) {
    z <- sis %>% filter(Sex == sx)
    data.frame(scheme = s, Sex = sx, cramers_V = cramer(z[[s]], z$Batch))
  }))
})) %>% pivot_wider(names_from = Sex, values_from = cramers_V)
print(ent %>% mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(),
      row.names = FALSE)

cat("\nSUS count per batch (stressed animals only):\n")
print(d %>% filter(Condition == "SIS") %>% group_by(Sex, Batch) %>%
        summarise(n = n(), yours = sum(yours == "SUS"), E2 = sum(E2 == "SUS"),
                  E3 = sum(E3 == "SUS"), .groups = "drop") %>% as.data.frame(),
      row.names = FALSE)

# ---------------------------------------------------------------------------
cat("\n======== 2. does batch actually shift the outcomes? ========\n")
cat("(if it does not, adding it only costs df; if it does, it buys precision)\n\n")
bx <- bind_rows(lapply(names(IND), function(nm) {
  m <- IND[[nm]]; z <- d %>% filter(!is.na(.data[[m]]))
  a <- anova(lm(reformulate("Batch", m), z))
  data.frame(outcome = nm, n = nrow(z),
             pct_var_batch = round(100 * a[1, "Sum Sq"] / sum(a[, "Sum Sq"]), 1),
             F = round(a[1, "F value"], 2), p = signif(a[1, "Pr(>F)"], 3))
}))
print(bx %>% arrange(desc(pct_var_batch)) %>% as.data.frame(), row.names = FALSE)

# ---------------------------------------------------------------------------
fit <- function(metric, labcol, groups, with_batch, stratum = NULL) {
  z <- d %>% filter(!is.na(.data[[metric]]), .data[[labcol]] %in% groups)
  if (!is.null(stratum)) z <- z %>% filter(Sex == stratum)
  z$g <- factor(z[[labcol]], levels = groups)
  if (dplyr::n_distinct(z$g) < 2) return(NULL)
  rhs <- if (with_batch && dplyr::n_distinct(z$Batch) > 1) c("g", "Batch") else "g"
  f <- lm(reformulate(rhs, metric), z)
  cn <- paste0("g", groups[2]); if (!cn %in% names(coef(f))) return(NULL)
  co <- summary(f)$coefficients[cn, ]; s <- summary(f)$sigma
  ci <- confint(f)[cn, ]
  # sigma  = residual SD. Adding Batch should SHRINK it -- that is the
  #          precision gain, and on the standardised scale it shows up as a
  #          larger d, not as a smaller SE(d): se_d = SE(coef)/sigma divides
  #          the gain straight back out. Report both so neither hides.
  # se_raw = SE of the raw coefficient. Its ratio mixes the sigma gain with
  #          the collinearity cost of a batch-entangled label.
  data.frame(outcome = names(IND)[match(metric, IND)], batch = with_batch,
             n = nrow(z), d = co[1] / s, lo = ci[1] / s, hi = ci[2] / s,
             se_d = co[2] / s, se_raw = co[2], sigma = s,
             t = co[3], p = co[4])
}

sweep <- function(labcol, groups, stratum = NULL) {
  raw <- bind_rows(lapply(unname(IND), function(m)
    bind_rows(fit(m, labcol, groups, FALSE, stratum),
              fit(m, labcol, groups, TRUE, stratum))))
  stopifnot(nrow(raw) > 0)
  raw %>%
    mutate(batch = ifelse(batch, "B", "noB")) %>%
    pivot_wider(names_from = batch,
                values_from = c(d, lo, hi, se_d, se_raw, sigma, t, p)) %>%
    mutate(sigma_ratio = sigma_B / sigma_noB,   # < 1 = batch removed noise
           seraw_ratio = se_raw_B / se_raw_noB, # < 1 = net precision gain
           t_ratio     = abs(t_B) / abs(t_noB)) %>%
    mutate(q_noB = p.adjust(p_noB, "BH"), q_B = p.adjust(p_B, "BH"))
}

shw <- function(x, ttl) {
  cat(sprintf("\n--- %s ---\n", ttl))
  print(x %>% mutate(d_noBatch = round(d_noB, 2), d_Batch = round(d_B, 2),
                     CI_Batch = sprintf("[%.2f, %.2f]", lo_B, hi_B),
                     sigma_r = round(sigma_ratio, 3),
                     seRaw_r = round(seraw_ratio, 3),
                     q_noBatch = signif(q_noB, 2), q_Batch = signif(q_B, 2)) %>%
          select(outcome, n, d_noBatch, d_Batch, CI_Batch, sigma_r, seRaw_r,
                 q_noBatch, q_Batch) %>% as.data.frame(), row.names = FALSE)
  cat(sprintf("  median: sigma ratio %.3f (%+.1f%% residual SD) | raw SE ratio %.3f (%+.1f%% SE)\n",
              median(x$sigma_ratio), 100 * (median(x$sigma_ratio) - 1),
              median(x$seraw_ratio), 100 * (median(x$seraw_ratio) - 1)))
}

cat("\n======== 3. CON vs SIS -- randomised, so batch is PURE precision ========\n")
cat("sigma_r / seRaw_r < 1 means adding Batch tightened the estimate.\n")
cs <- sweep("cond", c("CON", "SIS"))
shw(cs, "CON vs SIS (BH across the 9 independent outcomes)")
cat(sprintf("  median |change in d|: %.3f  -> the estimate itself is stable\n",
            median(abs(cs$d_B - cs$d_noB))))

cat("\n======== 4. SUS vs RES -- where the adjustment question bites ========\n")
for (s in c("yours", "E2", "E3")) {
  x <- sweep(s, c("RES", "SUS"))
  shw(x, sprintf("SUS vs RES, labels = %s", s))
  cat(sprintf("  median |d| without batch %.2f, with batch %.2f | median |change in d| %.3f\n",
              median(abs(x$d_noB)), median(abs(x$d_B)),
              median(abs(x$d_B - x$d_noB))))
}

cat("\n======== 5. sex-stratified (Batch = 3 levels inside each sex) ========\n")
for (sx in c("f", "m")) {
  for (s in c("yours", "E3")) {
    x <- sweep(s, c("RES", "SUS"), stratum = sx)
    shw(x, sprintf("SUS vs RES, labels = %s, %s only", s,
                   ifelse(sx == "f", "females", "males")))
  }
}

out <- bind_rows(lapply(c("yours", "E1", "E2", "E3"), function(s)
  sweep(s, c("RES", "SUS")) %>% mutate(scheme = s)))
write.csv(out, file.path(RES, "batch_in_model.csv"), row.names = FALSE)
cat("\nwrote:", file.path(RES, "batch_in_model.csv"), "\n")
