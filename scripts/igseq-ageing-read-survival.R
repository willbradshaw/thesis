###############################################################################
## FIGURE                                                                    ##
## Pilot study read survival plots                                           ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------
age_groups <- c("39", "56", "73", "128")

# Import counts table
tab <- import_counts(inpath)

# Extract replicate and individual identity
tab <- mutate(tab,
              INDIVIDUAL = sub("(\\d-\\d\\d)(.*)", "\\1", REPLICATE),
              REP = sub("(\\d-\\d\\d)(.*)", "\\2", REPLICATE),
              GROUP = as.numeric(sub("-.*", "", INDIVIDUAL)),
              AGE_DAYS = factor(age_groups[GROUP], levels = age_groups))

# Define palette
# TODO: Define age-group palette
palette <- c(colours_igseq[["ageing_group1"]], colours_igseq[["ageing_group2"]],
             colours_igseq[["ageing_group3"]], colours_igseq[["ageing_group4"]])
#------------------------------------------------------------------------------
# WRITE READ SURVIVAL COUNTS
#------------------------------------------------------------------------------

# Raw reads by replicate
nreads_raw_replicate <- tab %>% filter(STAGE == "presto_raw") %>%
  pull(NREADS_RAW)
nreads_raw_total <- round(sum(nreads_raw_replicate)/1e6, 1)
nreads_raw_replicate_min <- round(min(nreads_raw_replicate)/1e6, 1)
nreads_raw_replicate_max <- round(max(nreads_raw_replicate)/1e6, 1)

savetxt(nreads_raw_total, outpath_nreads_raw_total)
savetxt(nreads_raw_replicate_min, outpath_nreads_raw_replicate_min)
savetxt(nreads_raw_replicate_max, outpath_nreads_raw_replicate_max)

# Raw reads by individual
nreads_raw_individual <- tab %>% filter(STAGE == "presto_raw") %>%
  group_by(INDIVIDUAL) %>% summarise(NREADS_RAW = sum(NREADS_RAW)) %>%
  pull(NREADS_RAW)
nreads_raw_individual_min <- round(min(nreads_raw_individual)/1e6, 1)
nreads_raw_individual_max <- round(max(nreads_raw_individual)/1e6, 1)

savetxt(nreads_raw_individual_min, outpath_nreads_raw_individual_min)
savetxt(nreads_raw_individual_max, outpath_nreads_raw_individual_max)

# Percentage survival (makedb)
pcreads_surv_init <- tab %>% filter(STAGE == "changeo_makedb") %>% pull(NREADS_PC)
pcreads_surv_init_min <- round(min(pcreads_surv_init * 100), 1)
pcreads_surv_init_max <- round(max(pcreads_surv_init * 100), 1)
savetxt(pcreads_surv_init_min, outpath_surv_init_min)
savetxt(pcreads_surv_init_max, outpath_surv_init_max)

# Percentage lost during V-score filtering / clonotyping
pcreads_surv_all <- tab %>% filter(STAGE == "changeo_filter") %>% pull(NREADS_PC)
pcreads_surv_rel <- pcreads_surv_init - pcreads_surv_all
pcreads_surv_rel_min <- round(min(pcreads_surv_rel * 100), 1)
pcreads_surv_rel_max <- round(max(pcreads_surv_rel * 100), 1)
savetxt(pcreads_surv_rel_min, outpath_rel_loss_min)
savetxt(pcreads_surv_rel_max, outpath_rel_loss_max)

# Total lost during V-score filtering
total_lost_vscore <- tab %>% group_by(STAGE) %>% 
  summarise(NREADS = sum(NREADS)) %>% 
  filter(STAGE %in% c("changeo_makedb", "changeo_filter")) %>%
  pull(NREADS) %>% (function(x) abs(diff(x)))
pc_lost_vscore <- round(total_lost_vscore/sum(nreads_raw_replicate)*100, 1)
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
allplot_init <- countplot_all(tab, colour = "AGE_DAYS", stages_include = stages_init,
                              palette = palette, cname = "Age group (days)")
allplot_all <- countplot_all(tab, colour = "AGE_DAYS", stages_include = stages_all,
                              palette = palette, cname = "Age group (days)")

# Visualise plot
plot_unit = "cm"
plot_height<- 15
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
grid.draw(allplot_init)
popViewport(1)

# Save figure
plt <- grid.grab()
ggsave(plot = plt, filename = outpath_init, device = "svg", units = "cm",
       height=plot_height, width = plot_width)

vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
grid.draw(allplot_all)
popViewport(1)

# Save figure
plt <- grid.grab()
ggsave(plot = plt, filename = outpath_all, device = "svg", units = "cm",
       height=plot_height, width = plot_width)
