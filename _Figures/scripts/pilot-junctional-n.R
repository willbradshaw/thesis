###############################################################################
## FIGURE                                                                    ##
## Frequence of N characters in junctional sequences of pilot dataset        ##
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
tab_path <- "../_Data/changeo/ctabs/pilot_filtered.tab"

# Output paths
filename_nn <- "igseq-pilot-filtered-nn"
filename_1n <- "pilot-filtered-1n"

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(tab_path)

#------------------------------------------------------------------------------
# COUNT JUNCTIONAL N DISTRIBUTION
#------------------------------------------------------------------------------

# Separate counts for functional and nonfunctional
tab_nf <- tab %>% group_by(SEQUENCE_INPUT, JUNCTION, FUNCTIONAL) %>%
  summarise(CONSCOUNT = sum(CONSCOUNT)) %>% ungroup() %>%
  mutate(NN = str_count(JUNCTION, "N")) %>%
  group_by(FUNCTIONAL, NN) %>% 
  summarise(N = n(), CONSCOUNT = sum(CONSCOUNT)) %>% 
  group_by(FUNCTIONAL) %>% 
  mutate(N_PC = N/sum(N)*100, N_CUM = cumsum(N_PC), N_CUM = N_CUM-first(N_CUM))

# Collapse across functionality
tab_n <- tab_nf %>% group_by(NN) %>%
  summarise(N = sum(N), CONSCOUNT = sum(CONSCOUNT)) %>%
  mutate(N_PC = N/sum(N)*100, N_CUM = cumsum(N_PC), 
         N_CUM = N_CUM-first(N_CUM)) %>%
  filter(!is.na(NN)) %>%
  mutate(N_PC_N = N_PC/max(N_CUM) * 100, N_CUM_N = N_CUM/max(N_CUM) * 100)

# Convert into pretty table for printing
tab_n_pretty <- tab_n %>% select(`# junctional Ns` = NN,
                            `# unique sequences` = N,
                            `% of all sequences` = N_PC, 
                            `% of sequences with >0 junctional Ns` = N_PC_N)
tab_n_pretty[1,4] <- 0
tab_n_pretty[7,1] <- ">5"
tab_n_pretty[7,2] <- sum(tab_n_pretty[-(1:6), 2])
tab_n_pretty[7,3] <- sum(tab_n_pretty[-(1:6), 3])
tab_n_pretty[7,4] <- sum(tab_n_pretty[-(1:6), 4])
tab_n_pretty <- tab_n_pretty[1:7,]
tab_n_pretty[,3] <- signif(tab_n_pretty[,3], 3)
tab_n_pretty[,4] <- signif(tab_n_pretty[,4], 3)

# Save table
savetab(tab_n_pretty, filename_nn)

# Extract and save text values
nseq_nn_any <- tab_n %>% filter(NN > 0) %>% pull(N_PC) %>% sum() %>% floor()
nseq_1n_total <- tab_n %>% filter(NN==1) %>% pull(N_PC) %>% signif(3)
nseq_1n_withn <- tab_n %>% filter(NN==1) %>% pull(N_PC_N) %>% signif(3)
savetxt(nseq_nn_any, "pilot-filtered-nn-any")
savetxt(nseq_1n_total, paste0(filename_1n, "-total"))
savetxt(nseq_1n_withn, paste0(filename_1n, "-withn"))