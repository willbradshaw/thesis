###############################################################################
## FIGURE                                                                    ##
## Clone counts and statistics in gut dataset                                ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Parameters
age_groups = c("6", "16")
treatment_groups <- c("YI_6", "WT_16", "ABX_16", "SMT_16", "YMT_16")
palette <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_group2"]],
             colours_igseq[["gut_group3"]], colours_igseq[["gut_group4"]],
             colours_igseq[["gut_group5"]])
palette_age <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_allold"]])
palette_exp <- c(colours_igseq[["exp_pilot"]], colours_igseq[["exp_ageing"]],
                 colours_igseq[["exp_gut"]])
individuals_excluded <- c("1274", "1309")


#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab_gut <- import_tab(inpath_gut) %>% 
  mutate(GROUP = sub("YI_16", "YI_6", GROUP)) %>%
  filter(!INDIVIDUAL %in% individuals_excluded)
tab_age <- import_tab(inpath_age)
tab_pilot <- import_tab(inpath_pilot)

#------------------------------------------------------------------------------
# COUNT UNIQUE SEQUENCES WITH ASSIGNED CLONES
#------------------------------------------------------------------------------

nseq_naclone <- tab_gut %>% group_by(SEQUENCE_INPUT, CLONE) %>% summarise() %>% 
  pull(CLONE) %>% is.na %>% mean %>% (function(x) round((1-x)*100, 1))
savetxt(nseq_naclone, outpath_has_clone)

#------------------------------------------------------------------------------
# COMPARE CLONES (AND OTHER METRICS) PER INDIVIDUAL BETWEEN DATASETS
#------------------------------------------------------------------------------

# Get per-clone information
tab_cl_pilot <- tab_pilot %>% filter(!is.na(CLONE)) %>%
  group_by(AGE_DAYS, INDIVIDUAL, CLONE) %>%
  summarise(CLNCOUNT = n(), DUPCOUNT = sum(DUPCOUNT), 
            CONSCOUNT = sum(CONSCOUNT)) %>% mutate(EXPERIMENT = "pilot")

tab_cl_age <- tab_age %>% filter(!is.na(CLONE)) %>%
  group_by(AGE_DAYS, INDIVIDUAL, CLONE) %>%
  summarise(CLNCOUNT = n(), DUPCOUNT = sum(DUPCOUNT), 
            CONSCOUNT = sum(CONSCOUNT)) %>% mutate(EXPERIMENT = "ageing")

tab_cl_gut <- tab_gut %>% filter(!is.na(CLONE)) %>%
  group_by(AGE_WEEKS, GROUP, INDIVIDUAL, CLONE) %>%
  summarise(CLNCOUNT = n(), DUPCOUNT = sum(DUPCOUNT), 
            CONSCOUNT = sum(CONSCOUNT))  %>% mutate(EXPERIMENT = "gut") 

tab_cl_all <- bind_rows(tab_cl_pilot, tab_cl_age, tab_cl_gut) %>% 
  group_by(EXPERIMENT, INDIVIDUAL)

# Summarise into counts
tab_cl_all_counts <- tab_cl_all %>%
  summarise(N = n(), CLNCOUNT = sum(CLNCOUNT), DUPCOUNT = sum(DUPCOUNT), 
            CONSCOUNT = sum(CONSCOUNT))

# Combine, discard unshared variables, and melt 
tab_cl_all_counts_melted <- tab_cl_all_counts %>% ungroup %>%
  melt(id.vars = c("EXPERIMENT", "INDIVIDUAL")) %>%
  filter(variable != "CONSCOUNT") %>%
  mutate(EXPERIMENT = factor(EXPERIMENT, 
                             levels = c("pilot", "ageing", "gut")),
         variable = factor(variable,
                           levels = c("DUPCOUNT", "CLNCOUNT", "N"))
         )

# Make boxplots
theme_metrics <- 
  theme_classic() + theme_base + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(#size = fontsize_base * fontscale_title,
    angle = 45, hjust = 1),
    legend.title = element_text(margin = margin(r = 0.5, unit = "cm")),
    legend.text = element_text(margin = margin(r = 0.3, unit = "cm"))
  )
plot_metrics <- ggplot(tab_cl_all_counts_melted) + 
  geom_boxplot(aes(x=variable, y=value, fill = EXPERIMENT)) + 
  scale_x_discrete(breaks = c("DUPCOUNT", "CLNCOUNT", "N"),
                   labels = c("UMI groups", "Unique sequences", "Clones")) +
  scale_fill_manual(values = palette_exp, name = "Experiment") +
  ylab("Number in repertoire") + theme_metrics

plot_metrics_log <- plot_metrics + scale_y_log10()

# Extract and save text values for gut data
tab_cl_gut_counts <- tab_cl_all_counts %>% filter(EXPERIMENT == "gut")
clones_individual_min <- tab_cl_gut_counts %>% pull(N) %>% min
clones_individual_max <- tab_cl_gut_counts %>% pull(N) %>% max
clones_individual_med <- tab_cl_gut_counts %>% pull(N) %>% median %>% round
savetxt(clones_individual_min, outpath_clones_indiv_min)
savetxt(clones_individual_max, outpath_clones_indiv_max)
savetxt(clones_individual_med, outpath_clones_indiv_med)

# Test age and treatment effect on clone counts
tab_cl_counts_num <- tab_cl_gut %>% 
  group_by(AGE_WEEKS, GROUP, INDIVIDUAL, CLONE) %>%
  summarise(CLNCOUNT = sum(CLNCOUNT), DUPCOUNT = sum(DUPCOUNT), 
            CONSCOUNT = sum(CONSCOUNT)) %>%
  group_by(AGE_WEEKS, GROUP, INDIVIDUAL) %>% 
  summarise(N = n(), CLNCOUNT = sum(CLNCOUNT), DUPCOUNT = sum(DUPCOUNT), 
            CONSCOUNT = sum(CONSCOUNT)) %>% ungroup() %>%
  mutate(AGE_WEEKS = as.numeric(AGE_WEEKS))

cl_counts_kruskal_age <- kruskal.test(formula = N ~ AGE_WEEKS,
                                      data = tab_cl_counts_num)
cl_counts_kruskal_group <- kruskal.test(formula = N ~ as.factor(GROUP),
                                      data = tab_cl_counts_num)
savetxt(cl_counts_kruskal_age$p.value %>% signif(2), outpath_counts_kruskal_age)
savetxt(cl_counts_kruskal_group$p.value %>% signif(2), outpath_counts_kruskal_group)

# Plot range of clone numbers in different age/treatment groups
g_nclones_age <- ggplot(tab_cl_counts_num) + 
  geom_boxplot(aes(x=factor(AGE_WEEKS, levels = age_groups), y=N, 
                   fill = factor(AGE_WEEKS, levels = age_groups)),
               outlier.shape = NA) +
  geom_point(aes(x=factor(AGE_WEEKS, levels = age_groups), y=N), 
             size = 3, alpha = 0.5) +
  scale_fill_manual(values = palette_age, name = "Age group (weeks)") +
  xlab("Age at death (weeks)") +  scale_y_log10(name = "# Clones") +
  theme_classic() + theme_base + theme(legend.position = "none")
g_nclones_group <- ggplot(tab_cl_counts_num) + 
  geom_boxplot(aes(x=factor(GROUP, levels = treatment_groups), y=N, 
                   fill = factor(GROUP, levels = treatment_groups)),
               outlier.shape = NA) +
  geom_point(aes(x=factor(GROUP, levels = treatment_groups), y=N), 
             size = 3, alpha = 0.5) +
  scale_fill_manual(values = palette, name = "Treatment group") +
  scale_x_discrete(name = "Treatment group") +
  scale_y_log10(name = "# Clones") +
  theme_classic() + theme_base + theme(
    legend.position = "none",
    axis.text.x = element_text(size = fontsize_base * 0.8))

g_nclones <- gplot_grid(g_nclones_age, g_nclones_group)

ggsave(plot = g_nclones, filename = outpath_nclones, device = "svg", units = "cm",
       height=15, width = 25)

#------------------------------------------------------------------------------
# COMPUTE RELATIVE METRICS AND COMPARE BETWEEN EXPERIMENTS
#------------------------------------------------------------------------------

tab_cl_all_counts_relative <- tab_cl_all_counts %>%
  group_by(EXPERIMENT, INDIVIDUAL) %>%
  transmute(CLONES_PER_READ = N/CONSCOUNT, CLONES_PER_UMI = N/DUPCOUNT,
            CLONES_PER_SEQUENCE = N/CLNCOUNT)

tab_cl_all_counts_relative_melted <- tab_cl_all_counts_relative %>%
  melt(id.vars = c("EXPERIMENT", "INDIVIDUAL")) %>%
  filter(variable != "CLONES_PER_READ") %>%
  mutate(EXPERIMENT = factor(EXPERIMENT, 
                             levels = c("pilot", "ageing", "gut")),
         variable = factor(variable,
                           levels = c("CLONES_PER_UMI", "CLONES_PER_SEQUENCE"))
  )

plot_metrics_relative <- ggplot(tab_cl_all_counts_relative_melted,
                                aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill = EXPERIMENT), outlier.shape = NA) +
  scale_x_discrete(breaks = c("CLONES_PER_UMI", "CLONES_PER_SEQUENCE"),
                   labels = c("UMI groups", "Unique sequences")) +
  ylab("# per clone") +
  scale_fill_manual(values = palette_exp, name = "Experiment") + theme_metrics

g_metrics <- gplot_grid_onelegend(plot_metrics, plot_metrics_relative,
                                  plot_height = 15, plot_width = 25)
ggsave(plot = g_metrics, filename = outpath_comparative, device = "svg", units = "cm",
       height=15, width = 25)
