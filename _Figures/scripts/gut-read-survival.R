###############################################################################
## FIGURE                                                                    ##
## Gut study read survival plots                                             ##
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

# Configure input paths
survival_path <- "../_Data/changeo/survival/gut-reports-combined.tsv"
metadata_path <- "../_Data/fish_data/gut-samples.csv"

# Configure output
filename_raw <- "gut-reads-raw"
filename_base <- "gut-read-survival"
filename_init <- paste0(filename_base, "-init")
filename_all <- paste0(filename_base, "-all")
filename_rin <- paste0(filename_all, "-rin")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

# Import counts table
tab <- import_counts(survival_path) %>% mutate(INDIVIDUAL = REPLICATE) %>%
  full_join(suppressMessages(read_csv(metadata_path)), by = "INDIVIDUAL")

# Define palette
palette <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_group2"]],
             colours_igseq[["gut_group3"]], colours_igseq[["gut_group4"]],
             colours_igseq[["gut_group5"]])

#------------------------------------------------------------------------------
# WRITE READ SURVIVAL COUNTS
#------------------------------------------------------------------------------

# Raw reads
nreads_raw <- tab %>% filter(STAGE == "presto_raw") %>%
  pull(NREADS_RAW)
nreads_raw_total <- round(sum(nreads_raw)/1e6, 1)
nreads_raw_min <- round(min(nreads_raw)/1e6, 1)
nreads_raw_max <- round(max(nreads_raw)/1e6, 1)

savetxt(nreads_raw_total, paste0(filename_raw, "-total"))
savetxt(nreads_raw_min, paste0(filename_raw, "-min"))
savetxt(nreads_raw_max, paste0(filename_raw, "-max"))

# Percentage survival (makedb)
pcreads_surv_init <- tab %>% filter(STAGE == "changeo_makedb") %>% pull(NREADS_PC)
pcreads_surv_init_min <- round(min(pcreads_surv_init * 100), 1)
pcreads_surv_init_max <- round(max(pcreads_surv_init * 100), 1)
savetxt(pcreads_surv_init_min, paste0(filename_init, "-min"))
savetxt(pcreads_surv_init_max, paste0(filename_init, "-max"))

# Percentage lost during V-score filtering / clonotyping
pcreads_surv_all <- tab %>% filter(STAGE == "changeo_clonotyped") %>% 
  pull(NREADS_PC)
pcreads_surv_rel <- pcreads_surv_init - pcreads_surv_all
pcreads_surv_rel_min <- round(min(pcreads_surv_rel * 100), 1)
pcreads_surv_rel_max <- round(max(pcreads_surv_rel * 100), 1)
savetxt(pcreads_surv_rel_min, paste0(filename_base, "-rel-loss-min"))
savetxt(pcreads_surv_rel_max, paste0(filename_base, "-rel-loss-max"))

# Total lost during V-score filtering
total_lost_vscore <- tab %>% group_by(STAGE) %>% 
  summarise(NREADS = sum(NREADS)) %>% 
  filter(STAGE %in% c("changeo_makedb", "changeo_clonotyped")) %>%
  pull(NREADS) %>% (function(x) abs(diff(x)))
pc_lost_vscore <- round(total_lost_vscore/sum(nreads_raw)*100, 1)
savetxt(pc_lost_vscore, paste0(filename_base, "-rel-loss-total"))

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
savefig(plot = plt, filename = filename_init,
        height = plot_height, 
        width = plot_width)

vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
grid.draw(allplot_all)
popViewport(1)

# Save figure
plt <- grid.grab()
savefig(plot = plt, filename = filename_all,
        height = plot_height, 
        width = plot_width)

#------------------------------------------------------------------------------
# RE-PLOT SURVIVAL CURVES, COLOURED BY RIN
#------------------------------------------------------------------------------

# Specify colour palette
rin_colour_good <- "lightgray"
rin_colour_bad <- "black"

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

# Extract legend from absplot
g <- ggplotGrob(absplot_rin + theme(legend.position = "bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)

# Combine plots without legend
plt_rin <- plot_grid(absplot_rin + theme(legend.position = "none"),
                 relplot_rin + theme(legend.position = "none"), 
                 scaplot_rin + theme(legend.position = "none"),
                 ncol = 3, nrow = 1, labels="AUTO",
                 label_fontfamily = titlefont, label_fontface = "plain",
                 label_size = fontsize_base * fontscale_label)
combined_rin <- arrangeGrob(plt_rin, legend, ncol = 1, nrow = 2,
                            heights = unit.c(unit(1, "npc") - lheight, lheight))

# Draw and plot as before
plot_unit = "cm"
plot_height<- 15
plot_width <- 35
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
grid.draw(combined_rin)
popViewport(1)

# Save figure
plt_rin <- grid.grab()
savefig(plot = plt_rin, filename = filename_rin,
        height = plot_height, 
        width = plot_width)
