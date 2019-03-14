###############################################################################
## FIGURE                                                                    ##
## Depict ageing clones by abundance and overabundance                      ##
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
tab_path <- "../_Data/changeo/ctabs/ageing-final.tab"

# Output paths
filename_base <- "ageing-clones-expansions"

# Set parameters
palette <- c(colours_igseq[["ageing_group1"]], colours_igseq[["ageing_group2"]],
             colours_igseq[["ageing_group3"]], colours_igseq[["ageing_group4"]])
threshold_abs = 5
threshold_rel = 3
age_groups <- c("39", "56", "73", "128")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(tab_path)

#------------------------------------------------------------------------------
# COMPUTE CLONAL FREQUENCY AND OVERABUNDANCE
#------------------------------------------------------------------------------

# Basic clone table
clntab <- tab %>% group_by(INDIVIDUAL, SEQUENCE_INPUT, CLONE) %>%
  summarise() %>% compute_clntab() %>% arrange(desc(CLNFREQ)) %>%
  mutate(CLNFREQ = CLNFREQ * 100,
         CLNFREQ_REL = CLNFREQ / lead(CLNFREQ),
         CLNFREQ_REL = ifelse(is.na(CLNFREQ_REL), 1, CLNFREQ_REL))
clntab <- clntab %>% mutate(AGE_GROUP = as.integer(sub("-.*", "", INDIVIDUAL)),
                            AGE_DAYS = age_groups[AGE_GROUP])

#------------------------------------------------------------------------------
# PLOT AGAINST ROSENFELD THRESHOLDS
#------------------------------------------------------------------------------

g <- ggplot(clntab) +
  geom_point(aes(x=CLNFREQ, y=CLNFREQ_REL, colour = AGE_DAYS)) +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_x_continuous(name = "Relative clonal frequency (%)", 
                     limits = c(NA, max(clntab$CLNFREQ))) +
  scale_y_continuous(name = "Clonal overabundance (×)", 
                     limits = c(NA, max(clntab$CLNFREQ_REL))) +
  geom_vline(xintercept = threshold_abs, linetype = 2) +
  geom_hline(yintercept = threshold_rel, linetype = 2) +
  facet_wrap(~INDIVIDUAL, scales = "free") +
  theme_classic() + theme_base + theme(
    legend.position = "none",
    plot.margin = margin(t = 0.5, r = 0.5, l = 0.5, b = 0.5, unit = "cm")
  )

g_rep <- ggplot(clntab_rep) +
  geom_point(aes(x=CLNFREQ, y=CLNFREQ_REL, colour = INDIVIDUAL)) +
  facet_wrap(~REPLICATE, scales = "free", ncol = 3) +
  scale_colour_manual(values = palette, name = "Individual") +
  scale_x_continuous(name = "Relative clonal frequency (%)", 
                     limits = c(NA, max(clntab_rep$CLNFREQ))) +
  scale_y_continuous(name = "Clonal overabundance (×)", 
                     limits = c(NA, max(clntab_rep$CLNFREQ_REL))) +
  geom_vline(xintercept = threshold_abs, linetype = 2) +
  geom_hline(yintercept = threshold_rel, linetype = 2) +
  theme_classic() + theme_base + theme(
    legend.position = "none",
    plot.margin = margin(t = 0.5, r = 0.5, l = 0.5, b = 0.5, unit = "cm")
    )
  
#------------------------------------------------------------------------------
# SAVE PLOTS
#------------------------------------------------------------------------------

savefig(g, filename_base, height = 20, ratio = 1)
savefig(g_rep, paste0(filename_base, "-rep"), height = 20, ratio = 0.8)