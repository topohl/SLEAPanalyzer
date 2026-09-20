# ============================================================================
# 18_within_sex_composite.R
#
# THE PROBLEM
#
# The susceptibility composite is an equally-weighted average of six
# control-referenced z-scores. "Equally weighted" is only true if the six
# components have comparable spread. They do not, and the imbalance differs by
# sex. Variance share of the composite (DLSsingleSlim, SIS animals):
#
#             NOR   SucPref  dCORT  Bodyweight  Adrenal  Spleen
#   male     32.1%    15.3%  14.4%       18.2%    13.8%    6.3%
#   female    5.1%    12.0%  31.6%       21.0%    22.6%    7.6%
#
# So "susceptible" denotes a largely NOR-driven phenotype in males and a
# largely CORT/adrenal one in females. Same word, different construct -- and
# it is why the male NOR effect looked enormous (NOR is a third of the male
# composite) while the female one did not.
#
# THE FIX
#
# Centring is already sex-specific: components are centred on each batch's own
# controls and no batch is mixed-sex. Only the SCALING needs equalising. Each
# component is divided by its within-sex spread among the stressed animals, so
# every component contributes equally to the composite in both sexes.
#
# What this does and does not do:
#   - it makes the LABEL mean the same combination of symptoms in each sex
#   - it removes any real sex difference in how variable a given response is.
#     If females genuinely vary more in CORT response, that is discarded. That
#     is the deliberate trade: a consistent construct, not a faithful one.
#
# The control reference is kept for the centre, so a score still reads as
# "distance from this batch's controls" -- only the units change.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggplot2)
})

PROJ <- "C:/Users/topohl/iCloudDrive/Dokumente/Analysis/Behavior/correlate_sleap_boris"
SISD <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis"
ANA  <- "s:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis"
RES  <- file.path(PROJ, "results"); FIG <- file.path(PROJ, "figures")

canon <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(grepl("^[0-9]+$", x), sub("^0+(?=[0-9])", "", x, perl = TRUE), x)
}
SIX <- c("NOR", "SucPref", "dCORT", "Bodyweight", "Adrenal", "Spleen")

d <- suppressMessages(read_excel(file.path(SISD, "E9_Behavior_Data.xlsx"),
                                 sheet = "DLSsingleSlim")) %>%
  mutate(across(-any_of(c("ID", "Group", "Sex", "Batch")), ~suppressWarnings(as.numeric(.))),
         key = canon(ID), Batch = paste0("B", Batch))

sus_orig <- canon(trimws(readLines(file.path(ANA, "sus_animals.txt"), warn = FALSE)))

# --- equalise component weight within sex -----------------------------------
sis <- d %>% filter(Group == "SIS")
scale_by <- sis %>% group_by(Sex) %>%
  summarise(across(all_of(SIX), ~sd(., na.rm = TRUE)), .groups = "drop")
cat("within-sex SD used to rescale each component:\n")
print(scale_by %>% mutate(across(-Sex, ~round(., 2))) %>% as.data.frame(), row.names = FALSE)

eq <- sis
for (s in unique(sis$Sex)) {
  i <- eq$Sex == s
  for (cc in SIX) eq[[cc]][i] <- eq[[cc]][i] / scale_by[[cc]][scale_by$Sex == s]
}

sis$comp_current  <- rowMeans(sis[, SIX], na.rm = TRUE)
sis$comp_withinsex <- rowMeans(eq[, SIX], na.rm = TRUE)

vshare <- function(df, sx) {
  v <- sapply(SIX, function(cc) var(df[[cc]][df$Sex == sx], na.rm = TRUE))
  round(100 * v / sum(v), 1)
}
cat("\n=== variance share of the composite, before and after ===\n")
print(data.frame(component = SIX,
                 male_before = vshare(sis, "m"), male_after = vshare(eq, "m"),
                 female_before = vshare(sis, "f"), female_after = vshare(eq, "f")),
      row.names = FALSE)

# --- relabel ----------------------------------------------------------------
sis$L_current <- ifelse(sis$key %in% sus_orig, "SUS", "RES")
n_sus <- sum(sis$L_current == "SUS")
by_count <- function(z) ifelse(rank(z, ties.method = "first") <= n_sus, "SUS", "RES")
sis$L_withinsex <- by_count(sis$comp_withinsex)
sis$L_equalweight <- by_count(sis$comp_current)   # same threshold rule, old weighting

cat(sprintf("\n=== does the SUS assignment move? (holding the SUS count at %d) ===\n", n_sus))
cat(sprintf("  within-sex weighting vs your current list : %d of %d animals change\n",
            sum(sis$L_withinsex != sis$L_current), nrow(sis)))
print(table(current = sis$L_current, within_sex = sis$L_withinsex))
cat(sprintf("\n  SUS -> RES: %d   |   RES -> SUS: %d\n",
            sum(sis$L_current == "SUS" & sis$L_withinsex == "RES"),
            sum(sis$L_current == "RES" & sis$L_withinsex == "SUS")))

cat("\n=== changes by sex (is one sex relabelled more?) ===\n")
print(sis %>% group_by(Sex) %>%
        summarise(n = n(), SUS_current = sum(L_current == "SUS"),
                  SUS_withinsex = sum(L_withinsex == "SUS"),
                  changed = sum(L_current != L_withinsex), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\n=== animals that change ===\n")
print(sis %>% filter(L_current != L_withinsex) %>%
        mutate(move = paste(L_current, "->", L_withinsex),
               z_before = round(comp_current, 2), z_after = round(comp_withinsex, 2)) %>%
        select(ID, Batch, Sex, move, z_before, z_after) %>%
        arrange(Sex, move) %>% as.data.frame(), row.names = FALSE)

write.csv(sis %>% select(ID, key, Batch, Sex, comp_current, comp_withinsex,
                         L_current, L_withinsex),
          file.path(RES, "within_sex_composite_labels.csv"), row.names = FALSE)

# --- what it does to the CON/RES/SUS comparisons ---------------------------
dat <- read.delim(file.path(PROJ, "enriched", "sleap_all_batches_wide.tsv"),
                  stringsAsFactors = FALSE, na.strings = "") %>% mutate(key = canon(ID))
IND <- c("EPM open-arm fraction" = "sleap_EPM_open_frac",
         "EPM centre time"       = "sleap_EPM_center_s",
         "EPM distance"          = "sleap_EPM_distance_cm",
         "OFT centre %"          = "sleap_OFT_center_pct",
         "OFT corner %"          = "sleap_OFT_corner_pct",
         "SocP S1 preference"    = "sleap_SocP_S1_pref_index",
         "SocP S2 preference"    = "sleap_SocP_S2_pref_index")

j <- dat %>%
  left_join(sis %>% select(key, L_current, L_withinsex), by = "key") %>%
  mutate(P_current   = ifelse(Condition == "CON", "CON", L_current),
         P_withinsex = ifelse(Condition == "CON", "CON", L_withinsex))

contrast_sus_res <- function(lab, metric) {
  dd <- j %>% filter(!is.na(.data[[metric]]), .data[[lab]] %in% c("RES", "SUS"))
  if (nrow(dd) < 10) return(NULL)
  dd$g <- factor(dd[[lab]], levels = c("RES", "SUS"))
  f <- lm(reformulate(c("g", "Batch"), metric), dd)
  co <- summary(f)$coefficients["gSUS", ]
  data.frame(label = lab, metric = metric, d = co[1] / summary(f)$sigma, p = co[4])
}
cmp <- bind_rows(lapply(c("P_current", "P_withinsex"), function(l)
  bind_rows(lapply(unname(IND), function(m) contrast_sus_res(l, m))))) %>%
  mutate(label = recode(label, P_current = "current weighting",
                        P_withinsex = "within-sex equalised"),
         metric = names(IND)[match(metric, IND)]) %>%
  pivot_wider(names_from = label, values_from = c(d, p))

cat("\n=== SUS vs RES on the INDEPENDENT assays, under each weighting ===\n")
print(cmp %>% mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(),
      row.names = FALSE)
cat("\nThese assays are not in the classifier, so this is a fair comparison of\n",
    "whether the relabelling makes the groups behaviourally more distinct.\n", sep = "")
write.csv(cmp, file.path(RES, "within_sex_outcome_comparison.csv"), row.names = FALSE)

# --- Figure -----------------------------------------------------------------
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
pl <- bind_rows(
  data.frame(component = SIX, share = vshare(sis, "m"), Sex = "male", w = "current"),
  data.frame(component = SIX, share = vshare(eq, "m"),  Sex = "male", w = "within-sex equalised"),
  data.frame(component = SIX, share = vshare(sis, "f"), Sex = "female", w = "current"),
  data.frame(component = SIX, share = vshare(eq, "f"),  Sex = "female", w = "within-sex equalised")
) %>% mutate(component = factor(component, levels = rev(SIX)))
p <- ggplot(pl, aes(share, component, fill = Sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_vline(xintercept = 100 / 6, colour = AXIS, linetype = "22", linewidth = 0.5) +
  scale_fill_manual(values = c(male = "#2a78d6", female = "#eb6834"), name = NULL) +
  facet_wrap(~ w) +
  labs(title = "What each component contributes to the susceptibility score",
       subtitle = "Variance share within each sex. Dashed line is equal weight (1/6).",
       x = "share of composite variance (%)", y = NULL,
       caption = "Before equalising, NOR drives a third of the male score and a twentieth of the female one -- the same label meaning different things.") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = AXIS, linewidth = 0.4),
        axis.text = element_text(colour = MUTED, size = 8),
        axis.title = element_text(colour = INK2, size = 9),
        plot.title = element_text(colour = INK, face = "bold", size = 12),
        plot.subtitle = element_text(colour = INK2, size = 9),
        plot.caption = element_text(colour = MUTED, size = 8, hjust = 0),
        strip.text = element_text(colour = INK, face = "bold", size = 9),
        legend.position = "top")
ggsave(file.path(FIG, "fig12_within_sex_weighting.png"), p,
       width = 10, height = 4.4, dpi = 200, bg = SURFACE)

cat("\nwrote:", file.path(RES, "within_sex_composite_labels.csv"), "\n")
