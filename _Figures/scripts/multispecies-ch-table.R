###############################################################################
## FIGURE                                                                    ##
## Multispecies C-region maps                                                ##
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

# Configure input paths
tab_path <- "../_Data/species_data/multispecies_ch_regions.csv"

# Configure output
filename <- "multispecies-ch-regions"

#------------------------------------------------------------------------------
# READ TABLE
#------------------------------------------------------------------------------

tab <- suppressMessages(read_csv(tab_path)) %>%
  mutate(`\\textbf{Comments}` = gsub("&", "\\&", `\\textbf{Comments}`))

# More modifications?

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
  savetab(get(paste0("tab_", n)), paste0(filename, "-", n),
          align = align_string)
}