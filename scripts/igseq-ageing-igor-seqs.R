###############################################################################
## FIGURE                                                                    ##
## Pre-IGoR sequence counts in pilot dataset                                 ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

age_groups <- c("39", "56", "73", "128")

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

savetxt(min(counts), outpath_min)
savetxt(max(counts), outpath_max)

#------------------------------------------------------------------------------
# COUNT INDIVIDUAL SEQUENCES
#------------------------------------------------------------------------------

counts_group <- tab_unique %>% group_by(AGE_DAYS) %>% summarise(N = sum(N))
counts_group_out <- counts_group %>% 
  mutate(`Age group (days)` = factor(AGE_DAYS, levels = age_groups),
         `# Unique sequences` = N) %>% select(-AGE_DAYS, -N) %>%
  arrange(`Age group (days)`)

savetab(counts_group_out, outpath_tab)
