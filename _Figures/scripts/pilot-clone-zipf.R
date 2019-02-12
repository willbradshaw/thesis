###############################################################################
## FIGURE                                                                    ##
## FITTING RANK/ABUNDANCE DISTRIBUTIONS OF PILOT CLONES TO ZIPF'S LAW        ##
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

# Input paths
tab_path <- "../_Data/changeo/ctabs/pilot_clones.tab"

# Output paths
filename_base <- "pilot-clones-zipf"

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(tab_path)

#------------------------------------------------------------------------------
# COMPUTE CLONAL RANKS AND SIZES
#------------------------------------------------------------------------------

n_exclude <- c(6,3,2,5)
n_exclude_null <- c(0,0,0,0)
names(n_exclude) <- paste0("2-0", seq(3,6))
names(n_exclude_null) <- paste0("2-0", seq(3,6))

clntab <- tab %>% group_by(INDIVIDUAL, SEQUENCE_INPUT, CLONE) %>%
  summarise() %>%
  compute_clntab() %>% compute_expected_frequencies(n_exclude)

#------------------------------------------------------------------------------
# PLOT FIT AND RESIDUALS
#------------------------------------------------------------------------------

palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
             colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]])

zplot <- plot_zipf_fit(clntab) + 
  scale_colour_manual(values = palette, name = "Individual")
zplot_residual <- plot_zipf_residuals(clntab)

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

savefig(zplot, filename_base, height = 20, ratio = 1)