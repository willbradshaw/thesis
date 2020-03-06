###############################################################################
## FIGURE                                                                    ##
## Frequence of N characters in junctional sequences of pilot dataset        ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(inpath)

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
tab_n_pretty[,3] <- round(tab_n_pretty[,3], 3)
tab_n_pretty[,4] <- round(tab_n_pretty[,4], 2)

# Save table
savetab(tab_n_pretty, outpath_tab)

# Extract and save text values
nseq_nn_any <- tab_n %>% filter(NN > 0) %>% pull(N_PC) %>% sum() %>% floor()
nseq_1n_total <- tab_n %>% filter(NN==1) %>% pull(N_PC) %>% signif(3)
nseq_1n_withn <- tab_n %>% filter(NN==1) %>% pull(N_PC_N) %>% signif(3)
savetxt(nseq_nn_any, outpath_nn_any)
savetxt(nseq_1n_total, outpath_nseq_1n_total)
savetxt(nseq_1n_withn, outpath_nseq_1n_withn)
