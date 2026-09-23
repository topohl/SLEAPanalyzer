# ============================================================================
# 00_theme.R -- shared figure style for the Exp9 validation study
#
# Deliberately identical to `theme_exp9()` in
# SISanalyzer/exp9_publication/R/00_setup.R, so the validation figures and the
# publication figures read as one set: same palette, same type size, same
# gridline convention, same output widths.
#
# Nature conventions applied here:
#   - sans-serif, 7 pt body type
#   - 89 mm single-column / 120 mm / 183 mm double-column widths
#   - minimal non-data ink: no panel borders, no axis lines, no tick marks
#   - gridlines on the value axis only (override with grid = "both" for
#     scatters, "none" for heatmaps)
#   - vector PDF for submission plus a 600 dpi PNG preview from one call
#
# Palette: dark navy / neutral grey / coral. Navy and coral differ in hue AND
# lightness, so they stay separable under deuteranopia and protanopia (where
# the coral reads as tan), with grey sitting between them in lightness.
# ============================================================================

suppressMessages({ library(ggplot2); library(ggtext) })

MM  <- function(x) x / 25.4              # mm -> inches, for ggsave
W1  <- MM(89); W15 <- MM(120); W2 <- MM(183)

PAL <- c(CON = "#3F4576", RES = "#C2C2C2", SUS = "#F4636E",
         SIS = "#D98C93", Female = "#C0698F", Male = "#4A5490")

# Series colours for non-group figures (scatters, sweeps, Bland-Altman).
SER1 <- "#3F4576"   # navy  -- data
SER2 <- "#F4636E"   # coral -- fitted / emphasis
SER3 <- "#4A8C7A"   # teal  -- third series, when one is genuinely needed

SER4 <- "#C0698F"   # mauve -- fourth series

# Qualitative series palette for ordered variants (parameter sweeps, detector
# designs). Four hues that stay distinguishable in greyscale and under
# deuteranopia: navy, coral, teal, mauve -- drawn from the same family as the
# categorical group palette rather than a second colour language.
SERIES <- c(SER1, SER2, SER3, SER4)
INK  <- "#111111"; MUTED <- "#767676"; AXIS <- "#5E5E5E"; GRID <- "#EAEAEA"
RULE <- "#CFCFCF"

# Diverging scale for correlation matrices: navy (negative) to coral
# (positive) through near-white, so it shares the categorical palette's hues
# instead of introducing a second colour language.
DIV_NEG <- "#3F4576"; DIV_MID <- "#F4F4F2"; DIV_POS <- "#F4636E"

darken <- function(hex, f = 0.55) {
  rgb <- grDevices::col2rgb(hex) * f
  grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
}
PAL_DARK <- setNames(darken(PAL), names(PAL))

theme_exp9 <- function(base = 7, grid = c("y", "x", "both", "none")) {
  grid <- match.arg(grid)
  gl <- element_line(colour = GRID, linewidth = 0.3)
  theme_minimal(base_size = base, base_family = "sans") +
    theme(
      text               = element_text(colour = INK),
      axis.text          = element_text(colour = MUTED, size = base - 0.5),
      axis.title         = element_text(colour = AXIS, size = base),
      axis.ticks         = element_blank(),
      axis.line          = element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "both")) gl else element_blank(),
      panel.grid.major.x = if (grid %in% c("x", "both")) gl else element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_blank(),
      panel.spacing      = unit(7, "pt"),
      strip.background   = element_blank(),
      strip.text         = element_markdown(colour = INK, size = base, hjust = 0,
                                            face = "bold", margin = margin(b = 2.5)),
      legend.position    = "top",
      legend.key.size    = unit(6, "pt"),
      legend.text        = element_text(size = base - 0.5),
      legend.title       = element_text(size = base - 0.5, colour = AXIS),
      legend.margin      = margin(0, 0, 0, 0),
      legend.background  = element_blank(),
      plot.title         = element_text(size = base + 1.5, face = "bold", hjust = 0),
      plot.subtitle      = element_markdown(size = base - 0.5, colour = AXIS,
                                            hjust = 0, lineheight = 1.4),
      plot.caption       = element_markdown(size = base - 1.5, colour = MUTED,
                                            hjust = 0, lineheight = 1.4),
      plot.tag           = element_text(size = base + 2, face = "bold"),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA)
    )
}

# Coloured group names instead of a keyed legend, as in the publication set.
# A visible separator, not repeated spaces: markdown collapses whitespace runs.
colour_key <- function(groups = c("CON", "RES", "SUS"))
  paste(sprintf("<span style='color:%s'>**%s**</span>", PAL_DARK[groups], groups),
        collapse = " <span style='color:#C8C8C8'>|</span> ")

# Vector for the journal, raster for quick viewing, from one call.
save_fig <- function(plot, name, width, height, dir = FIG) {
  ggsave(file.path(dir, paste0(name, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 600, bg = "white")
  invisible(name)
}
