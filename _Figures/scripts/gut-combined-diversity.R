###############################################################################
## FIGURE                                                                    ##
## Plot diversity spectra of gut clonal repertoires                          ##
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

# Configure input
# TODO: Select input settings to match other spectra in chapter
seqset <- "all" # or "functional"
copy <- "NULL" # or "DUPCOUNT"
groupby <- "group" # or "age"

tab_path_grouped_clone <- "../_Data/changeo/spectra/gut-age_clone-diversity-grouped_seqs-all_copy-NULL.tsv"
tab_path_grouped_vj <- "../_Data/changeo/spectra/gut-age_VJ-diversity-grouped_seqs-all_copy-NULL.tsv"
tab_path_solo_clone <- "../_Data/changeo/spectra/gut-age_clone-diversity-solo_seqs-all_copy-NULL.tsv"
tab_path_solo_vj <- "../_Data/changeo/spectra/gut-age_VJ-diversity-solo_seqs-all_copy-NULL.tsv"

# Set parameters
palette <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_allold"]])
age_groups <- c("6", "16")
qvals <- c(0,1,1.5,2,3,4)
significance_level <- 0.05

# Output paths
filename_base <- "igseq-gut-combined-diversity"

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab_grouped_clone <- import_div(tab_path_grouped_clone) %>% mutate(TYPE = "Clonal repertoire")
tab_grouped_vj <- import_div(tab_path_grouped_vj) %>% mutate(TYPE = "VJ repertoire")
tab_grouped <- bind_rows(tab_grouped_clone, tab_grouped_vj)

tab_solo_clone <- import_div(tab_path_solo_clone) %>% mutate(TYPE = "Clonal repertoire")
tab_solo_vj <- import_div(tab_path_solo_vj) %>% mutate(TYPE = "VJ repertoire")
tab_solo <- bind_rows(tab_solo_clone, tab_solo_vj)

# Process grouped spectra for alpha plot
tab_grouped <- mutate(tab_grouped, 
                      AGE_WEEKS = factor(sub(".*_", "", AGE_WEEKS),
                                         levels = age_groups))

# Process solo spectra as appropriate
tab_solo <- mutate(tab_solo,
                         AGE_WEEKS = factor(sub(".*_", "", AGE_WEEKS),
                                            levels = age_groups))

#------------------------------------------------------------------------------
# GENERATE AND SAVE ALPHA SPECTRA
#------------------------------------------------------------------------------

g_alpha <- plot_diversity_alpha(tab_grouped, "AGE_WEEKS") +
  scale_colour_manual(values = palette, name = "Age at death (weeks)") +
  scale_fill_manual(values = palette, name = "Age at death (weeks)") +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm"))) +
  facet_wrap(.~TYPE, scales = "free_y")

savefig(plot = g_alpha, filename = paste0(filename_base, "-alpha"),
        height = 15, width = 25)