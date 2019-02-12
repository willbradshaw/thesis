###############################################################################
## FIGURE                                                                    ##
## Plot diversity spectra of pilot clonal repertoires                        ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

rm(list = ls())

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")
source("aux/ggplot2.R")
source("aux/changeo.R")

# Configure input
# TODO: Select input settings to match other spectra in chapter
seqset <- "all" # or "functional"
copy <- "NULL" # or "DUPCOUNT"
tab_path <- paste0("../_Data/changeo/spectra/pilot_clone-diversity-grouped_",
                   "seqs-", seqset, "_copy-", copy, ".tsv")

# Output paths
filename_base <- "pilot-clone-diversity"

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_div(tab_path)

#------------------------------------------------------------------------------
# GENERATE ALPHA AND BETA SPECTRA
#------------------------------------------------------------------------------

palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
             colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]])

g_alpha <- plot_diversity_alpha(tab, "INDIVIDUAL") +
  scale_colour_manual(values = palette, name = "Individual") +
  scale_fill_manual(values = palette, name = "Individual")
g_beta <- plot_diversity_beta_scaled(tab, "INDIVIDUAL") +
  scale_colour_manual(values = palette, name = "Individual") +
  scale_fill_manual(values = palette, name = "Individual")

#------------------------------------------------------------------------------
# COMBINE SPECTRA WITH SINGLE LEGEND
#------------------------------------------------------------------------------

# Extract legend from absplot
g <- ggplotGrob(g_alpha)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)
# Combine plots without legend
plt <- plot_grid(g_alpha + theme(legend.position = "none"),
                 g_beta + theme(legend.position = "none"), 
                 ncol = 2, nrow = 1, labels="AUTO",
                 label_fontfamily = titlefont, label_fontface = "plain",
                 label_size = fontsize_base * fontscale_label)
combined <- arrangeGrob(plt,
                        legend,
                        ncol = 1, nrow = 2,
                        heights = unit.c(unit(1, "npc") - lheight, lheight))

# Visualise plot
plot_unit = "cm"
plot_height <- 15
plot_width <- 25
map_layout <- grid.layout(
  ncol = 1,
  nrow = 1,
  heights = unit(plot_height, plot_unit),
  widths = unit(plot_width, plot_unit)
)
vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
grid.draw(combined)
popViewport(1)

plt <- grid.grab()

#------------------------------------------------------------------------------
# SAVE FIGURES
#------------------------------------------------------------------------------

savefig(plot = plt, filename = filename_base,
        height = plot_height, 
        width = plot_width)