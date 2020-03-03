###############################################################################
## FIGURE                                                                    ##
## Gut study read survival plots                                             ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Parameters
treatment_groups <- c("YI_6", "WT_16", "ABX_16", "SMT_16", "YMT_16")
palette <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_group2"]],
             colours_igseq[["gut_group3"]], colours_igseq[["gut_group4"]],
             colours_igseq[["gut_group5"]])
rin_colour_good <- "blue"
rin_colour_bad <- "red"

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

# Import counts table
tab <- import_counts(inpath_survival) %>% mutate(INDIVIDUAL = REPLICATE) %>%
  full_join(suppressMessages(read_csv(inpath_metadata)), by = "INDIVIDUAL") %>%
  mutate(GROUP = factor(sub("wk", "", GROUP), levels = treatment_groups))

#------------------------------------------------------------------------------
# WRITE READ SURVIVAL COUNTS
#------------------------------------------------------------------------------

# Raw reads
nreads_raw <- tab %>% filter(STAGE == "presto_raw") %>%
  pull(NREADS_RAW)
nreads_raw_total <- round(sum(nreads_raw)/1e6, 1)
nreads_raw_min <- round(min(nreads_raw)/1e6, 1)
nreads_raw_max <- round(max(nreads_raw)/1e6, 1)

savetxt(nreads_raw_total, outpath_nreads_raw_total)
savetxt(nreads_raw_min, outpath_nreads_raw_min)
savetxt(nreads_raw_max, outpath_nreads_raw_max)

# Percentage survival (makedb)
pcreads_surv_init <- tab %>% filter(STAGE == "changeo_makedb") %>% pull(NREADS_PC)
pcreads_surv_init_min <- round(min(pcreads_surv_init * 100), 1)
pcreads_surv_init_max <- round(max(pcreads_surv_init * 100), 1)
savetxt(pcreads_surv_init_min, outpath_surv_init_min)
savetxt(pcreads_surv_init_max, outpath_surv_init_max)

# Percentage lost during V-score filtering / clonotyping
pcreads_surv_all <- tab %>% filter(STAGE == "changeo_clonotyped") %>% 
  pull(NREADS_PC)
pcreads_surv_rel <- pcreads_surv_init - pcreads_surv_all
pcreads_surv_rel_min <- round(min(pcreads_surv_rel * 100), 1)
pcreads_surv_rel_max <- round(max(pcreads_surv_rel * 100), 1)
savetxt(pcreads_surv_rel_min, outpath_rel_loss_min)
savetxt(pcreads_surv_rel_max, outpath_rel_loss_max)

# Total lost during V-score filtering
total_lost_vscore <- tab %>% group_by(STAGE) %>% 
  summarise(NREADS = sum(NREADS)) %>% 
  filter(STAGE %in% c("changeo_makedb", "changeo_clonotyped")) %>%
  pull(NREADS) %>% (function(x) abs(diff(x)))
pc_lost_vscore <- round(total_lost_vscore/sum(nreads_raw)*100, 1)
savetxt(pc_lost_vscore, outpath_rel_loss_total)

#------------------------------------------------------------------------------
# PLOT INITIAL (UP TO CHANGEO_MAKE) SURVIVAL CURVES
#------------------------------------------------------------------------------

# Define stages
stages_init <- c("presto_raw", "presto_filter", "presto_mask", "presto_correct",
                 "presto_merged", "presto_collapsed", "presto_split",
                 "changeo_makedb")
stages_all <- c(stages_init, "changeo_filter", "changeo_clonotyped")

# Make grid table for plot
allplot_init <- countplot_all(tab, colour = "GROUP", stages_include = stages_init,
                              palette = palette, cname = "Treatment group")
allplot_all <- countplot_all(tab, colour = "GROUP", stages_include = stages_all,
                             palette = palette, cname = "Treatment group",
                             hline_rel = 0.3)
ggsave(plot = allplot_init, filename = outpath_init, device = "svg", units = "cm",
       height=15, width = 25)
ggsave(plot = allplot_all, filename = outpath_all, device = "svg", units = "cm",
       height=15, width = 25)

#------------------------------------------------------------------------------
# RE-PLOT SURVIVAL CURVES, COLOURED BY RIN
#------------------------------------------------------------------------------

# Compute correlation between RIN and read survival
tab_rin <- tab %>% filter(STAGE == "changeo_clonotyped") %>% 
  select(INDIVIDUAL, RIN, LP, NREADS_PC) %>% 
  mutate(PC_READS_LOST = pcreads_surv_rel)
ct <- cor.test(~RIN + NREADS_PC, data = tab_rin)
cor_label <- paste0("r = ", signif(ct$estimate, 3), 
                "\np = ", signif(ct$p.value, 3))

# Make simple RIN/read survival scatter plot
scaplot_rin <- ggplot(tab_rin) +
  geom_point(aes(x=RIN, y=NREADS_PC, colour = RIN), size = 3) +
  annotate("text", label = cor_label, x = 5.8, y=0.05, family = font,
           size = 4) +
  xlab("RNA integrity number") + 
  ylab("% Read survival") +
  scale_y_continuous(labels = function(y) y * 100, limits = c(0,1)) +
  scale_color_continuous(low = rin_colour_bad, high = rin_colour_good,
                         name = "RNA integrity number") +
  theme_classic() + theme_base

# Make absolute and relative plots
absplot_rin <-   countplot_base("RIN") + countline_abs(tab) + 
  scale_y_continuous(name = "# Reads", limits = c(0, max(tab$NREADS))) +
  countscale_x(stages_include = stages_all) +
  scale_color_continuous(low = rin_colour_bad, high = rin_colour_good,
                         name = "RNA integrity number")
relplot_rin <-   countplot_base("RIN") + countline_rel(tab) + 
  scale_y_continuous(name = "% Reads", limits = c(0, max(tab$NREADS))) +
  countscale_x(stages_include = stages_all) +
  scale_color_continuous(low = rin_colour_bad, high = rin_colour_good,
                         name = "RNA integrity number") +
  geom_hline(yintercept = 0.3, colour = "red", linetype = 2)

# Make and save plot
plt_rin <- gplot_grid_onelegend(absplot_rin, relplot_rin, scaplot_rin,
                                nrow = 1, ncol = 3, plot_height = 15,
                                plot_width = 35)

ggsave(plot = plt_rin, filename = outpath_rin, device = "svg", units = "cm",
       height=15, width = 35)
