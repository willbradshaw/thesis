###############################################################################
## FIGURE                                                                    ##
## Combined nfu/xma RSS composition plots                                    ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")
source("aux/ggplot2.R")

# Source single-species plots
source("scripts/nfu-rss-seqlogo.R")
source("scripts/xma-rss-seqlogo.R")

# Configure output
filename_combined <- "combined-rss-seqlogo-all"

#------------------------------------------------------------------------------
# MAKE COMBINED PLOTS
#------------------------------------------------------------------------------

plt_combined <- plot_grid(ggHeptamer_nfu$all + ggtitle("Heptamer composition"),
             spacer_all_annotated_nfu + ggtitle("Spacer length"),
             ggNonamer_nfu$all + ggtitle("Nonamer composition"),
             ggHeptamer_xma$all + ggtitle("Heptamer composition"),
             spacer_all_annotated_xma + ggtitle("Spacer length"),
             ggNonamer_xma$all + ggtitle("Nonamer composition"),
             ncol = 3, nrow = 2, labels="AUTO", label_fontfamily = titlefont,
             label_fontface = "plain",
             label_size = fontsize_base * fontscale_label)

#------------------------------------------------------------------------------
# SAVE PLOTS
#------------------------------------------------------------------------------

plot_height <- 20
plot_ratio <- 3/2 # Width vs height

savefig(plot = plt_combined, filename = filename_combined, device = "svg", 
        height = plot_height, ratio = plot_ratio)
savefig(plot = plt_combined, filename = filename_combined
        , device = "png", 
        height = plot_height, ratio = plot_ratio)