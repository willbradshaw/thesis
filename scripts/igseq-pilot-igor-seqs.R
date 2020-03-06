###############################################################################
## FIGURE                                                                    ##
## Pre-IGoR sequence counts in pilot dataset                                 ##
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
# COUNT SEQUENCES
#------------------------------------------------------------------------------

counts <- tab %>% group_by(INDIVIDUAL) %>% summarise(N = n()) %>% pull(N)

savetxt(min(counts), outpath_min)
savetxt(max(counts), outpath_max)
savetxt(sum(counts), outpath_sum)
