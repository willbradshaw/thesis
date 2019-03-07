###############################################################################
## FIGURES & DATA                                                            ##
## Pilot study repertoire composition                                        ##
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

# Input path
tab_path <- "../_Data/changeo/ctabs/pilot-named.tab"
# TODO: Change to a path in the thesis directory

# Configure output
filename_nseq <- "pilot-nseq-init"
filename_nreads <- "pilot-nreads-init"
filename_vmean <- "pilot-vscore-mean"
filename_vsd <- "pilot-vscore-sd"
filename_base <- "pilot-functional"
filename_proportion <- paste0(filename_base, "-prop")
filename_vscores <- paste0(filename_base, "-vscores")

#------------------------------------------------------------------------------
# COUNT UNIQUE SEQUENCES
#------------------------------------------------------------------------------

tab <- import_tab(tab_path)

# Number of unique sequences per replicate
n_replicate <- tab %>% group_by(REPLICATE) %>% count

# Number of unique sequences per individual
n_individual <- tab %>% group_by(INDIVIDUAL, SEQUENCE_INPUT) %>% count %>%
  group_by(INDIVIDUAL) %>% count

# Total number of unique sequences
n_total <- tab %>% group_by(SEQUENCE_INPUT) %>% count %>% 
  ungroup() %>% count

# Extract values to save
nseq_total <- n_total %>% pull(nn)
nseq_individual_min <- n_individual %>% 
  pull(nn) %>% (function(x) floor(min(x)/10^3)*10^3)
nseq_individual_max <- n_individual %>% 
  pull(nn) %>% (function(x) floor(max(x)/10^3)*10^3)
nseq_replicate_min <- n_replicate %>% 
  pull(n) %>% (function(x) floor(min(x)/10^3)*10^3)
nseq_replicate_max <- n_replicate %>% 
  pull(n) %>% (function(x) floor(max(x)/10^3)*10^3)

# Save values
savetxt(nseq_total, paste0(filename_nseq, "-total"))
savetxt(nseq_individual_min, paste0(filename_nseq, "-individual-min"))
savetxt(nseq_individual_max, paste0(filename_nseq, "-individual-max"))
savetxt(nseq_replicate_min, paste0(filename_nseq, "-replicate-min"))
savetxt(nseq_replicate_max, paste0(filename_nseq, "-replicate-max"))

#------------------------------------------------------------------------------
# ANALYSE FUNCTIONAL COMPOSITION
#------------------------------------------------------------------------------

# Diagnose functionality of each sequence
tab_func <- diagnose_functionality(tab %>% mutate(HAS_J = !is.na(J_CALL)))

# Count functional classes in pooled repertoire
tab_countfunc <- tab_func %>% group_by(FUNC_STATE, FUNC_DESC, SEQUENCE_INPUT,
                                       FUNCTIONAL, HAS_J) %>% 
  summarise(CONSCOUNT = sum(CONSCOUNT), DUPCOUNT = sum(DUPCOUNT), 
            V_SCORE = mean(V_SCORE)) # TODO: Include more columns?

#------------------------------------------------------------------------------
# DESCRIBE FUNCTIONAL COMPOSITION
#------------------------------------------------------------------------------

tab_countfunc_summ <- tab_countfunc %>% group_by(FUNCTIONAL, HAS_J) %>%
  summarise(N = n(), CONSCOUNT = sum(CONSCOUNT), DUPCOUNT = sum(DUPCOUNT),
            VSCORE_MEAN = mean(V_SCORE), VSCORE_SD = sd(V_SCORE)) %>%
  ungroup() %>% 
  mutate(N_PC = N/sum(N)*100,
         CONSCOUNT_PC = CONSCOUNT/sum(CONSCOUNT)*100,
         DUPCOUNT_PC = DUPCOUNT/sum(DUPCOUNT)*100)

# Percentage of unique sequences
nseq_functional <- tab_countfunc_summ %>% filter(FUNCTIONAL) %>% 
  pull(N_PC) %>% round(1)
nseq_noj <- tab_countfunc_summ %>% filter(!HAS_J) %>% 
  pull(N_PC) %>% round(1)
nseq_other <- tab_countfunc_summ %>% filter(!FUNCTIONAL, HAS_J) %>%
  pull(N_PC) %>% round(1)
savetxt(nseq_functional, paste0(filename_nseq, "-pc-functional"))
savetxt(nseq_noj, paste0(filename_nseq, "-pc-noj"))
savetxt(nseq_other, paste0(filename_nseq, "-pc-other"))

# Percentage of reads
nreads_functional <- tab_countfunc_summ %>% filter(FUNCTIONAL) %>% 
  pull(CONSCOUNT_PC) %>% round(1)
nreads_noj <- tab_countfunc_summ %>% filter(!HAS_J) %>%
  pull(CONSCOUNT_PC) %>% round(1)
nreads_other <- tab_countfunc_summ %>% filter(!FUNCTIONAL, HAS_J) %>%
  pull(CONSCOUNT_PC) %>% round(1)
savetxt(nreads_functional, paste0(filename_nreads, "-pc-functional"))
savetxt(nreads_noj, paste0(filename_nreads, "-pc-noj"))
savetxt(nreads_other, paste0(filename_nreads, "-pc-other"))

# Mean V-scores
vmean_functional <- tab_countfunc_summ %>% filter(FUNCTIONAL) %>% 
  pull(VSCORE_MEAN) %>% signif(3)
vmean_noj <- tab_countfunc_summ %>% filter(!HAS_J) %>%
  pull(VSCORE_MEAN) %>% signif(3)
vmean_other <- tab_countfunc_summ %>% filter(!FUNCTIONAL, HAS_J) %>%
  pull(VSCORE_MEAN) %>% signif(3)
savetxt(vmean_functional, paste0(filename_vmean, "-functional"))
savetxt(vmean_noj, paste0(filename_vmean, "-noj"))
savetxt(vmean_other, paste0(filename_vmean, "-other"))

# V-score SDs
vsd_functional <- tab_countfunc_summ %>% filter(FUNCTIONAL) %>% 
  pull(VSCORE_SD) %>% signif(3)
vsd_noj <- tab_countfunc_summ %>% filter(!HAS_J) %>%
  pull(VSCORE_SD) %>% signif(3)
vsd_other <- tab_countfunc_summ %>% filter(!FUNCTIONAL, HAS_J) %>%
  pull(VSCORE_SD) %>% signif(3)
savetxt(vsd_functional, paste0(filename_vsd, "-functional"))
savetxt(vsd_noj, paste0(filename_vsd, "-noj"))
savetxt(vsd_other, paste0(filename_vsd, "-other"))

#------------------------------------------------------------------------------
# PLOT FUNCTIONAL COMPOSITION OF UNFILTERED SEQUENCES
#------------------------------------------------------------------------------

tab_countfunc_melt <- tab_countfunc %>% 
  group_by(FUNC_STATE, FUNC_DESC) %>%
  summarise(N = n(), CONSCOUNT = sum(CONSCOUNT), DUPCOUNT = sum(DUPCOUNT)) %>%
  ungroup() %>%
  melt(variable.name = c("VAR"), value.name = "VALUE", 
       id.var = c("FUNC_STATE", "FUNC_DESC")) %>% 
  group_by(VAR) %>%
  mutate(PC = VALUE / sum(VALUE) * 100)

tab_countfunc_vscore <- tab_func %>% group_by(FUNCTIONAL, HAS_J) %>% 
  summarise(VSCORE_MEAN = mean(V_SCORE), VSCORE_SD = sd(V_SCORE))

func_levels <- tab_countfunc_melt %>% arrange(FUNC_STATE) %>% pull(FUNC_DESC) %>% unique

plot_proportions <- ggplot(tab_countfunc_melt) + 
  geom_col(aes(x=factor(VAR, levels = c("CONSCOUNT", "DUPCOUNT", "N")), 
               y=PC, fill = factor(FUNC_DESC, levels = func_levels)),
           position = "stack") + 
  scale_x_discrete(breaks = c("CONSCOUNT", "DUPCOUNT", "N"), 
                   labels = c("Reads", "UMI groups", "Unique sequences")) + 
  scale_fill_discrete(name = "State") + ylab("%") + 
  theme_base + theme(axis.title.x = element_blank(), 
                     axis.text.x = element_text(size = fontsize_base * fontscale_title,
                                                angle = 45, hjust = 1),
                     legend.justification = "center")

plot_vscores <- ggplot(tab_countfunc) +
  geom_density(aes(x=V_SCORE, 
                   colour = factor(FUNC_DESC, levels = func_levels))) + 
  facet_grid(FUNC_DESC~., scales = "free_y") + 
  scale_colour_discrete(name = "State") + xlab("V-alignment score") + ylab("Density") +
  theme_base + theme(
    axis.text.y = element_blank(), strip.text.y = element_blank()
  )

#------------------------------------------------------------------------------
# FILTER BY V SCORE AND RE-ANALYSE
#------------------------------------------------------------------------------

# Diagnose functionality of each sequence
tab_func_filtered <- tab_func %>%
  filter(V_SCORE >= 100)

# Count functional classes in pooled repertoire
tab_countfunc_filtered <- tab_func_filtered %>% 
  group_by(FUNC_STATE, FUNC_DESC, SEQUENCE_INPUT, FUNCTIONAL, HAS_J) %>% 
  summarise(CONSCOUNT = sum(CONSCOUNT), DUPCOUNT = sum(DUPCOUNT), 
            V_SCORE = mean(V_SCORE)) # TODO: Include more columns?

tab_countfunc_filtered_melt <- tab_countfunc_filtered %>%
  group_by(FUNC_STATE, FUNC_DESC) %>%
  summarise(N = n(), CONSCOUNT = sum(CONSCOUNT), DUPCOUNT = sum(DUPCOUNT)) %>%
  ungroup() %>%
  melt(variable.name = c("VAR"), value.name = "VALUE", 
       id.var = c("FUNC_STATE", "FUNC_DESC")) %>% as.tibble() %>%
  group_by(VAR) %>%
  mutate(PC = VALUE / sum(VALUE) * 100)

plot_proportions_filtered <- ggplot(tab_countfunc_filtered_melt) + 
  geom_col(aes(x=factor(VAR, levels = c("CONSCOUNT", "DUPCOUNT", "N")), 
               y=PC, fill = factor(FUNC_DESC, levels = func_levels))) + 
  scale_x_discrete(breaks = c("CONSCOUNT", "DUPCOUNT", "N"), 
                   labels = c("Reads", "UMI groups", "Unique sequences")) + 
  scale_fill_discrete(name = "State") + ylab("%") + 
  theme_base + theme(axis.title.x = element_blank(), 
                     axis.text.x = element_text(size = fontsize_base * fontscale_title,
                                                angle = 45, hjust = 1),
                     legend.justification = "center")


# Number and percentage of reads post-filtering
tab_countfunc_filtered_summ <- tab_countfunc_filtered %>% group_by(FUNCTIONAL, HAS_J) %>%
  summarise(N = n(), CONSCOUNT = sum(CONSCOUNT), DUPCOUNT = sum(DUPCOUNT),
            VSCORE_MEAN = mean(V_SCORE), VSCORE_SD = sd(V_SCORE)) %>%
  ungroup() %>% 
  mutate(N_PC = N/sum(N)*100,
         CONSCOUNT_PC = CONSCOUNT/sum(CONSCOUNT)*100,
         DUPCOUNT_PC = DUPCOUNT/sum(DUPCOUNT)*100)

nseq_filtered_functional <- tab_countfunc_filtered_summ %>% filter(FUNCTIONAL) %>% 
  pull(N) %>% round(1)
nseq_filtered_noj <- tab_countfunc_filtered_summ %>% filter(!HAS_J) %>%
  pull(N) %>% round(1)
nseq_filtered_other <- tab_countfunc_filtered_summ %>% filter(!FUNCTIONAL, HAS_J) %>%
  pull(N) %>% round(1)
nseq_filtered_functional_pc <- tab_countfunc_filtered_summ %>% filter(FUNCTIONAL) %>% 
  pull(N_PC) %>% round(1)
nseq_filtered_noj_pc <- tab_countfunc_filtered_summ %>% filter(!HAS_J) %>%
  pull(N_PC) %>% round(1)
nseq_filtered_other_pc <- tab_countfunc_filtered_summ %>% filter(!FUNCTIONAL, HAS_J) %>%
  pull(N_PC) %>% round(1)
savetxt(nseq_filtered_functional, paste0(filename_nseq, "-functional-filtered"))
savetxt(nseq_filtered_noj, paste0(filename_nseq, "-noj-filtered"))
savetxt(nseq_filtered_other, paste0(filename_nseq, "-other-filtered"))
savetxt(nseq_filtered_functional_pc, paste0(filename_nseq, "-pc-functional-filtered"))
savetxt(nseq_filtered_noj_pc, paste0(filename_nseq, "-pc-noj-filtered"))
savetxt(nseq_filtered_other_pc, paste0(filename_nseq, "-pc-other-filtered"))

#------------------------------------------------------------------------------
# COUNT SEQUENCES REMOVED IN THIS WAY
#------------------------------------------------------------------------------

nseq_dropped <- nrow(tab_countfunc) - nrow(tab_countfunc_filtered)
nreads_dropped <- tab_countfunc_filtered %>% pull(CONSCOUNT) %>% sum %>%
  (function(x) x/(tab_countfunc %>% pull(CONSCOUNT) %>% sum)) %>%
  (function(y) round((1-y)*100, 1))
nseq_kept <- nrow(tab_countfunc_filtered)
nreads_kept <- tab_countfunc_filtered %>% pull(CONSCOUNT) %>% sum %>%
  (function(x) x/(tab_countfunc %>% pull(CONSCOUNT) %>% sum)) %>%
  (function(y) round(y*100, 1))

savetxt(nseq_kept, paste0(filename_nseq, "-kept-vscore"))
savetxt(nreads_kept, paste0(filename_nreads, "-pc-kept-vscore"))
savetxt(nseq_dropped, paste0(filename_nseq, "-dropped-vscore"))
savetxt(nreads_dropped, paste0(filename_nreads, "-pc-dropped-vscore"))

#------------------------------------------------------------------------------
# COMBINE PROPORTION PLOTS WITH SINGLE LEGEND
#------------------------------------------------------------------------------

# Extract legend from absplot
g <- ggplotGrob(plot_proportions + theme(legend.position = "bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)
# Combine plots without legend
plt <- plot_grid(plot_proportions + ggtitle("V-score ≥ 0") + theme(legend.position = "none"),
                plot_proportions_filtered + ggtitle("V-score ≥ 100") + theme(legend.position = "none"), 
                ncol = 2, nrow = 1, labels="AUTO",
                label_fontfamily = titlefont, label_fontface = "plain",
                label_size = fontsize_base * fontscale_label)
combined <- arrangeGrob(plt,
                        legend,
                        ncol = 1, nrow = 2,
                        heights = unit.c(unit(1, "npc") - lheight, lheight))

# Visualise plot
plot_unit = "cm"
plot_height <- 15
plot_width <- 25
map_layout <- grid.layout(
  ncol = 1,
  nrow = 1,
  heights = unit(plot_height, plot_unit),
  widths = unit(plot_width, plot_unit)
)
vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
grid.draw(combined)
popViewport(1)

plt <- grid.grab()

#------------------------------------------------------------------------------
# SAVE FIGURES
#------------------------------------------------------------------------------

savefig(plot = plt, filename = filename_proportion,
        height = plot_height, 
        width = plot_width)
savefig(plot = plot_vscores, filename = filename_vscores,
        height = plot_height, width = plot_width)
