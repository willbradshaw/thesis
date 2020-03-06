###############################################################################
## FIGURE                                                                    ##
## Multispecies C-region maps                                                ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# READ TABLE
#------------------------------------------------------------------------------

tab <- suppressMessages(read_csv(inpath)) %>%
  mutate(`\\textbf{Comments}` = gsub("&", "\\&", `\\textbf{Comments}`))

#------------------------------------------------------------------------------
# WRITE TABLES
#------------------------------------------------------------------------------

align_string <- paste0(">{\\italic}", 
                       paste(rep("l", ncol(tab)),collapse = ""))

# Split table to fit on page
rows_per_tab <- 25
tab_split_at <- seq(0, ceiling(nrow(tab)/rows_per_tab)*rows_per_tab, 
                    rows_per_tab)
row_n_tab <- length(tab_split_at) - 1
for (n in 1:row_n_tab){
  assign(paste0("tab_", n), tab[(tab_split_at[n]+1):tab_split_at[n+1],])
  savetab(get(paste0("tab_", n)), get(paste0("outpath_ch", n)),
          align = align_string)
}
