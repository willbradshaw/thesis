###############################################################################
## FIGURE                                                                    ##
## Clone counts and statistics in gut dataset                                ##
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
inpath_pilot <- "../_Data/changeo/ctabs/pilot-final.tab"
inpath_age <- "../_Data/changeo/ctabs/ageing-final.tab"
inpath_gut <- "../_Data/changeo/ctabs/gut-final.tab"

# Output paths
filename_comparative <- "igseq-comparative"
filename_nclones <- "igseq-gut-nclones"
filename_base <- "igseq-gut-clone-sizes"

# Parameters
palette <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_group2"]],
             colours_igseq[["gut_group3"]], colours_igseq[["gut_group4"]],
             colours_igseq[["gut_group5"]])
treatment_groups <- c("ABX_16", "SMT_16", "WT_16", "YI_6", "YMT_16")
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
# COUNT UNIQUE SEQUENCES WITH NA CLONES
#------------------------------------------------------------------------------

nseq_naclone <- tab_gut %>% group_by(SEQUENCE_INPUT, CLONE) %>% summarise() %>% 
  pull(CLONE) %>% is.na %>% mean %>% (function(x) round((1-x)*100, 1))
#savetxt(nseq_naclone, "gut-nseq-assigned-clones")

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
savetxt(clones_individual_min, paste0(filename_nclones, "-individual-min"))
savetxt(clones_individual_max, paste0(filename_nclones, "-individual-max"))
savetxt(clones_individual_med, paste0(filename_nclones, "-individual-med"))

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

plot_metrics_relative <- ggplot(tab_cl_all_counts_relative_melted) + 
  geom_boxplot(aes(x=variable, y=value, fill = EXPERIMENT)) + 
  scale_x_discrete(breaks = c("CLONES_PER_UMI", "CLONES_PER_SEQUENCE"),
                   labels = c("UMI groups", "Unique sequences")) +
  ylab(expression(Clones ^ -1 )) +
  scale_fill_manual(values = palette_exp, name = "Experiment") + theme_metrics

#------------------------------------------------------------------------------
# SAVE FIGURES
#------------------------------------------------------------------------------

# Save metric boxplots with single shared legend

# Extract legend
g <- ggplotGrob(plot_metrics)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)
# Combine plots without legend
plt_metrics <- plot_grid(plot_metrics + theme(legend.position = "none"),
                         plot_metrics_relative + theme(legend.position = "none"), 
                         ncol = 2, nrow = 1, labels="AUTO",
                         label_fontfamily = titlefont, label_fontface = "plain",
                         label_size = fontsize_base * fontscale_label)
metrics_combined <- arrangeGrob(plt_metrics, legend, ncol = 1, nrow = 2,
                        heights = unit.c(unit(1, "npc") - lheight, lheight))

# Visualise plot
plot_height = 15
plot_width = 25
plot_unit = "cm"
map_layout <- grid.layout(ncol = 1, nrow = 1,
  heights = unit(plot_height, plot_unit), widths = unit(plot_width, plot_unit)
)
vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
grid.draw(metrics_combined)
popViewport(1)

plt_metrics <- grid.grab()
savefig(plt_metrics, filename = paste0(filename_comparative, "-metrics"),
        height = plot_height, width = plot_width)

