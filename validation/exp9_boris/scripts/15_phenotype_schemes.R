# ============================================================================
# 15_phenotype_schemes.R
#
# Two questions.
#
# A. DOES THE SUS LIST CHANGE under different labelling schemes, and how would
#    we know which one to trust? There is no ground truth for susceptibility,
#    so "which label is right" cannot be answered directly. What CAN be
#    measured is STABILITY: each batch z-score is referenced to only 4 control
#    animals, so the label of any animal near the threshold depends on which
#    4 controls it happened to be compared against. Bootstrapping the control
#    reference turns that into a per-animal flip rate.
#
# B. DOES IT MATTER for the outcome? Crossing the label scheme with whether
#    batch is in the outcome model separates the two corrections, which fix
#    different problems: z-scoring fixes misclassification of the PREDICTOR,
#    batch-as-covariate fixes noise in the OUTCOME.
#
# Classifier components (Analysis/SIS_Analysis/E9_Behavior_Data.xlsx,
# DLSsingleSlim): NOR D2, sucrose preference, dCORT, body weight, adrenal,
# spleen. NOR is an input, so every NOR outcome here uses a LEAVE-NOR-OUT
# label -- see 14_leave_one_out_phenotype.R.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggplot2)
})

PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
SIS  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis"
ANA  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis"
RES  <- file.path(PROJ, "results"); FIG <- file.path(PROJ, "figures")

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

# --- z-scoring machinery ----------------------------------------------------
# Centre on each batch's OWN control mean -- that is the batch effect being
# removed. Scale either by that batch's control SD (n = 4, unstable) or by the
# pooled within-sex control SD (n = 12, the same logic with a firmer estimate).
zs <- function(d, col, ctrl_idx, pooled_sd) {
  out <- rep(NA_real_, nrow(d))
  for (b in unique(d$Batch)) {
    ib <- which(d$Batch == b)
    ci <- intersect(ib, ctrl_idx)
    cv <- d[[col]][ci]; cv <- cv[is.finite(cv)]
    if (length(cv) < 2) next
    s <- if (pooled_sd) {
      sx <- d$Sex[ib][1]
      pv <- d[[col]][intersect(which(d$Sex == sx), ctrl_idx)]
      sd(pv[is.finite(pv)])
    } else sd(cv)
    if (!is.finite(s) || s == 0) next
    out[ib] <- (d[[col]][ib] - mean(cv)) / s
  }
  out
}
composite <- function(d, cols, ctrl_idx, pooled_sd) {
  rowMeans(vapply(cols, function(cc) zs(d, cc, ctrl_idx, pooled_sd),
                  numeric(nrow(d))), na.rm = TRUE)
}

ctrl_idx <- which(comp$Group == "CON")
sis_idx  <- which(comp$Group == "SIS")

comp$z6_batch  <- composite(comp, SIX,  ctrl_idx, FALSE)
comp$z6_pooled <- composite(comp, SIX,  ctrl_idx, TRUE)
comp$z5_batch  <- composite(comp, FIVE, ctrl_idx, FALSE)
comp$z5_pooled <- composite(comp, FIVE, ctrl_idx, TRUE)

sis <- comp[sis_idx, ] %>% filter(is.finite(z6_batch), is.finite(z5_pooled))
sis$L_orig <- ifelse(sis$key %in% sus_orig, "SUS", "RES")
sis$L_bc   <- ifelse(sis$key %in% sus_bc,   "SUS", "RES")
n_sus <- sum(sis$L_orig == "SUS")

# The existing batch-corrected list is reproduced almost perfectly by a cutoff
# on z6_batch, so recover that cutoff rather than inventing one. A FIXED cutoff
# (not a quantile) is used from here: a quantile would force the same SUS rate
# in every batch and erase real between-batch differences in susceptibility.
cand <- sort(unique(sis$z6_batch))
acc <- vapply(cand, function(k) mean((sis$z6_batch <= k) == (sis$L_bc == "SUS")), numeric(1))
CUT <- cand[which.max(acc)]
cat(sprintf("implied z cutoff reproducing sus_animals_batchCorrected: %.3f (accuracy %.1f%%)\n",
            CUT, 100 * max(acc)))

lab <- function(z) ifelse(z <= CUT, "SUS", "RES")
sis$L_z_batch   <- lab(sis$z6_batch)
sis$L_z_pooled  <- lab(sis$z6_pooled)
sis$L_loo_batch <- lab(sis$z5_batch)
sis$L_loo_pool  <- lab(sis$z5_pooled)

SCHEMES <- c(L_orig = "uncorrected roster", L_bc = "batchCorrected (existing)",
             L_z_batch = "z-scored, per-batch SD", L_z_pooled = "z-scored, pooled SD",
             L_loo_batch = "leave-NOR-out, per-batch SD",
             L_loo_pool = "leave-NOR-out, pooled SD")

# --- A. do the labels change? ----------------------------------------------
cat(sprintf("\nSIS animals: %d | SUS under the uncorrected roster: %d\n", nrow(sis), n_sus))
cat("\n=== A1. how many animals are SUS under each scheme, and how many differ from the roster ===\n")
tab <- bind_rows(lapply(names(SCHEMES), function(s)
  data.frame(scheme = SCHEMES[[s]], n_SUS = sum(sis[[s]] == "SUS"),
             differs_from_roster = sum(sis[[s]] != sis$L_orig))))
print(as.data.frame(tab), row.names = FALSE)

cat("\n=== A2. SUS rate per batch (a fixed threshold on a shifting scale shows up here) ===\n")
print(sis %>% group_by(Batch, Sex) %>%
        summarise(n = n(), across(all_of(names(SCHEMES)), ~sum(. == "SUS")), .groups = "drop") %>%
        rename_with(~ifelse(.x %in% names(SCHEMES), SCHEMES[.x], .x)) %>%
        as.data.frame(), row.names = FALSE)

# --- A3. stability: bootstrap the control reference -------------------------
# The honest way to ask "would this animal be labelled the same way again".
# Resample the CONTROL animals within each batch, rebuild the composite, and
# relabel. An animal that flips often is one the threshold cannot resolve.
set.seed(1)
B <- 400
boot_flip <- function(cols, pooled_sd) {
  base <- lab(composite(comp, cols, ctrl_idx, pooled_sd))[sis_idx][
    match(sis$key, comp$key[sis_idx])]
  flips <- matrix(NA, nrow = nrow(sis), ncol = B)
  for (b in seq_len(B)) {
    ci <- unlist(lapply(split(ctrl_idx, comp$Batch[ctrl_idx]),
                        function(v) sample(v, length(v), replace = TRUE)))
    z <- composite(comp, cols, ci, pooled_sd)
    l <- lab(z)[match(sis$key, comp$key)]
    flips[, b] <- l != base
  }
  rowMeans(flips, na.rm = TRUE)
}
sis$flip_batch  <- boot_flip(SIX, FALSE)
sis$flip_pooled <- boot_flip(SIX, TRUE)

cat(sprintf("\n=== A3. label stability under %d bootstraps of the control reference ===\n", B))
cat(sprintf("  per-batch SD (n=4 controls) : mean flip rate %.1f%% | animals flipping >10%% of the time: %d of %d\n",
            100 * mean(sis$flip_batch), sum(sis$flip_batch > 0.10), nrow(sis)))
cat(sprintf("  pooled within-sex SD (n=12) : mean flip rate %.1f%% | animals flipping >10%% of the time: %d of %d\n",
            100 * mean(sis$flip_pooled), sum(sis$flip_pooled > 0.10), nrow(sis)))

cat("\n--- the least stable animals (per-batch SD) ---\n")
print(sis %>% arrange(desc(flip_batch)) %>% head(10) %>%
        mutate(z6 = round(z6_batch, 2), flip = sprintf("%.0f%%", 100 * flip_batch)) %>%
        select(ID, Batch, Sex, z6, roster = L_orig, z_label = L_z_batch, flip) %>%
        as.data.frame(), row.names = FALSE)

write.csv(sis %>% select(ID, key, Batch, Sex, all_of(names(SCHEMES)),
                         z6_batch, z6_pooled, z5_batch, z5_pooled,
                         flip_batch, flip_pooled),
          file.path(RES, "phenotype_scheme_labels.csv"), row.names = FALSE)

# --- B. does it matter for the outcome? -------------------------------------
dat <- read.delim(file.path(PROJ, "enriched", "sleap_all_batches_wide.tsv"),
                  stringsAsFactors = FALSE, na.strings = "") %>% mutate(key = canon(ID))
j <- dat %>% inner_join(sis %>% select(key, all_of(names(SCHEMES))), by = "key")

fit_one <- function(label_col, metric, use_batch) {
  d <- j %>% filter(!is.na(.data[[metric]]), .data[[label_col]] %in% c("RES", "SUS"))
  if (nrow(d) < 10) return(NULL)
  d$g <- factor(d[[label_col]], levels = c("RES", "SUS"))
  fit <- lm(reformulate(if (use_batch) c("g", "Batch") else "g", metric), data = d)
  co <- summary(fit)$coefficients["gSUS", ]; ci <- confint(fit)["gSUS", ]
  s <- summary(fit)$sigma
  data.frame(scheme = SCHEMES[[label_col]], batch_in_model = use_batch, metric = metric,
             n = nrow(d), d = co[1] / s, lo = ci[1] / s, hi = ci[2] / s, p = co[4],
             stringsAsFactors = FALSE)
}

# NOR outcomes use leave-NOR-out labels only. The NOR-containing schemes are
# shown for NON-NOR outcomes, where they are not circular.
nor_rows <- bind_rows(lapply(c("L_loo_batch", "L_loo_pool"), function(l)
  bind_rows(lapply(c(TRUE, FALSE), function(ub) fit_one(l, "sleap_NOR_D2", ub)))))
oth_rows <- bind_rows(lapply(names(SCHEMES), function(l)
  bind_rows(lapply(c("sleap_SocP_S1_pref_index", "sleap_EPM_open_frac", "sleap_OFT_center_pct"),
                   function(m) bind_rows(lapply(c(TRUE, FALSE), function(ub) fit_one(l, m, ub)))))))
out <- bind_rows(nor_rows, oth_rows)
write.csv(out, file.path(RES, "phenotype_scheme_outcomes.csv"), row.names = FALSE)

cat("\n=== B. NOR D2 (leave-NOR-out labels only): z-scoring x batch-in-model ===\n")
print(nor_rows %>% mutate(d = round(d, 2), CI = sprintf("[%.2f, %.2f]", lo, hi),
                          p = signif(p, 3)) %>%
        select(scheme, batch_in_model, n, d, CI, p) %>% as.data.frame(), row.names = FALSE)

cat("\n=== B2. assays NOT in the classifier (never circular) ===\n")
print(oth_rows %>% filter(batch_in_model, scheme %in% SCHEMES[c("L_orig", "L_z_pooled")]) %>%
        mutate(d = round(d, 2), p = signif(p, 3)) %>%
        select(scheme, metric, n, d, p) %>% as.data.frame(), row.names = FALSE)

# --- Figure: stability ------------------------------------------------------
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
pl <- sis %>%
  select(z6_batch, `per-batch SD (n=4)` = flip_batch, `pooled within-sex SD (n=12)` = flip_pooled) %>%
  pivot_longer(-z6_batch, names_to = "reference", values_to = "flip")
p <- ggplot(pl, aes(z6_batch, flip, colour = reference)) +
  geom_vline(xintercept = CUT, colour = AXIS, linetype = "22", linewidth = 0.5) +
  geom_point(size = 1.6, alpha = 0.8) +
  scale_colour_manual(values = c("#eb6834", "#2a78d6"), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "How reliable is each animal's SUS/RES label?",
       subtitle = sprintf("%d bootstraps of the control reference. Dashed line is the classification threshold.", B),
       x = "susceptibility composite (batch z-scored)",
       y = "how often the label flips",
       caption = "Animals near the threshold flip most. Widening the control reference from 4 to 12 animals is what reduces it.") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = AXIS, linewidth = 0.4),
        axis.text = element_text(colour = MUTED, size = 8),
        axis.title = element_text(colour = INK2, size = 9),
        plot.title = element_text(colour = INK, face = "bold", size = 12),
        plot.subtitle = element_text(colour = INK2, size = 9),
        plot.caption = element_text(colour = MUTED, size = 8, hjust = 0),
        legend.position = "top")
ggsave(file.path(FIG, "fig10_label_stability.png"), p,
       width = 9, height = 4.6, dpi = 200, bg = SURFACE)

cat("\nwrote:", file.path(RES, "phenotype_scheme_labels.csv"), "\n")
