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

tab_unique <- tab %>% group_by(AGE_DAYS, INDIVIDUAL, SEQUENCE_INPUT) %>%
  summarise(N = 1)
counts <- tab_unique %>% group_by(INDIVIDUAL) %>% summarise(N = sum(N)) %>%
  pull(N)

savetxt(min(counts), paste0(filename_base, "-min"))
savetxt(max(counts), paste0(filename_base, "-max"))

#------------------------------------------------------------------------------
# COUNT INDIVIDUAL SEQUENCES
#------------------------------------------------------------------------------

counts_group <- tab_unique %>% group_by(AGE_DAYS) %>% summarise(N = sum(N))
counts_group_out <- counts_group %>% 
  mutate(`Age group (days)` = factor(AGE_DAYS, levels = age_groups),
         `# Unique sequences` = N) %>% select(-AGE_DAYS, -N) %>%
  arrange(`Age group (days)`)

savetab(counts_group_out, filename_base)