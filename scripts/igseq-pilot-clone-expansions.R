###############################################################################
## FIGURE                                                                    ##
## Depict pilot clones by abundance and overabundance                        ##
###############################################################################
aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Set parameters
palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
             colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]])
threshold_abs = 5
threshold_rel = 3

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(inpath)

#------------------------------------------------------------------------------
# COMPUTE CLONAL FREQUENCY AND OVERABUNDANCE
#------------------------------------------------------------------------------

# Basic clone table
clntab <- tab %>% group_by(INDIVIDUAL, SEQUENCE_INPUT, CLONE) %>%
  summarise() %>% compute_clntab() %>% arrange(desc(CLNFREQ)) %>%
  mutate(CLNFREQ = CLNFREQ * 100,
         CLNFREQ_REL = CLNFREQ / lead(CLNFREQ),
         CLNFREQ_REL = ifelse(is.na(CLNFREQ_REL), 1, CLNFREQ_REL))

clntab_rep <- tab %>% group_by(INDIVIDUAL=REPLICATE, SEQUENCE_INPUT, CLONE) %>%
  summarise() %>% compute_clntab() %>% arrange(desc(CLNFREQ)) %>%
  mutate(CLNFREQ = CLNFREQ * 100,
         CLNFREQ_REL = CLNFREQ / lead(CLNFREQ),
         CLNFREQ_REL = ifelse(is.na(CLNFREQ_REL), 1, CLNFREQ_REL),
         REPLICATE = sub("(2-0\\d)", "\\1 ", INDIVIDUAL)) %>% 
  ungroup() %>% mutate(INDIVIDUAL = sub(" .*", "", REPLICATE)) %>% 
  group_by(INDIVIDUAL)

#------------------------------------------------------------------------------
# PLOT AGAINST ROSENFELT THRESHOLDS
#------------------------------------------------------------------------------

g <- ggplot(clntab) +
  geom_point(aes(x=CLNFREQ, y=CLNFREQ_REL, colour = INDIVIDUAL)) +
  facet_wrap(~INDIVIDUAL, scales = "free") +
  scale_colour_manual(values = palette, name = "Individual") +
  scale_x_continuous(name = "Relative clonal frequency (%)", 
                     limits = c(NA, max(clntab$CLNFREQ))) +
  scale_y_continuous(name = "Clonal overabundance (×)", 
                     limits = c(NA, max(clntab$CLNFREQ_REL))) +
  geom_vline(xintercept = threshold_abs, linetype = 2) +
  geom_hline(yintercept = threshold_rel, linetype = 2) +
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

ggsave(plot = g, filename = outpath_all, device = "svg", units = "cm",
       height=20, width = 20)
ggsave(plot = g_rep, filename = outpath_rep, device = "svg", units = "cm",
       height=20, width = 20*0.8)
