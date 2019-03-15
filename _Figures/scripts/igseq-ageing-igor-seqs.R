###############################################################################
## FIGURE                                                                    ##
## Pre-IGoR sequence counts in pilot dataset                                 ##
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
source("aux/changeo.R")

# Specify input path
experiment <- "ageing"
inpath <- paste0("../_Data/changeo/ctabs/", experiment, "-igor.tab")

# Specify output path
filename_base <- paste0("igseq-", experiment, "-igor-seqs")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(inpath)

#------------------------------------------------------------------------------
# COUNT INDIVIDUAL SEQUENCES
#------------------------------------------------------------------------------

counts <- tab %>% group_by(INDIVIDUAL) %>% summarise(N = n()) %>% pull(N)

savetxt(min(counts), paste0(filename_base, "-min"))
savetxt(max(counts), paste0(filename_base, "-max"))
savetxt(sum(counts), paste0(filename_base, "-sum"))

#------------------------------------------------------------------------------
# COUNT INDIVIDUAL SEQUENCES
#------------------------------------------------------------------------------

counts_group <- tab %>% group_by(AGE_DAYS) %>% summarise(N = n())

# savetxt(min(counts), paste0(filename_base, "-min"))
# savetxt(max(counts), paste0(filename_base, "-max"))
# savetxt(sum(counts), paste0(filename_base, "-sum"))