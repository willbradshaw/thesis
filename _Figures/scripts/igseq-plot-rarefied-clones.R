###############################################################################
## ANALYSIS                                                                  ##
## Rarefaction plots of clonal counts for different datasets                 ##
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
inpath <- "../_Data/changeo/rarefied_clones.tab"

# Specify parameters
scale <- "DUPCOUNT" # or NULL
clone_field <- "CLONE"
palette_exp <- c(colours_igseq[["exp_pilot"]], colours_igseq[["exp_ageing"]],
                 colours_igseq[["exp_gut"]])
P <- 20
plot_height <- 15
plot_width <- 25

# Output paths
outpath <- "igseq-rarefied-clone"

#------------------------------------------------------------------------------
# READ INPUT
#------------------------------------------------------------------------------

col <- cols(MEAN="d", SD="d", SAMPLE_SIZE="d", .default="c")
r <- import_tsv(inpath, col = col) %>%
  mutate(EXPERIMENT = factor(EXPERIMENT, levels = c("pilot", "age", "gut")))

#------------------------------------------------------------------------------
# PLOT RAREFACTION CURVES OF CLONE COUNTS
#------------------------------------------------------------------------------

# Define axis titles
x_unit <- list("NULL" = "unique sequences",
               "DUPCOUNT" = "UMI groups",
               "CONSCOUNT" = "input reads")
x_title <- paste("Number of", x_unit[[scale]])
y_unit <- list("CLONE" = "clones") # TODO: Expand as needed
y_title <- paste("Number of", y_unit[[clone_field]])

# Plot mean lines only
plot_nclones_means <- ggplot(r %>% filter(METRIC == "N_CLONES")) +
  geom_line(aes(x=SAMPLE_SIZE, y=MEAN, colour = EXPERIMENT, 
                group = interaction(EXPERIMENT,INDIVIDUAL))) +
  scale_color_manual(values = palette_exp, name = "Experiment") +
  scale_x_continuous(breaks = seq(0, 10000, 2000), 
                     labels = function(x) x/1000, 
                     name = paste(x_title, "(×1000)")) +
  scale_y_continuous(labels = function(y) y/1000, 
                     name = paste(y_title, "(×1000)")) +
  theme_classic() + theme_base

# Add SD regions
plot_nclones_sd <- plot_nclones_means +
  geom_ribbon(aes(x=SAMPLE_SIZE, ymin=MEAN-SD, ymax = MEAN + SD, 
                  fill = EXPERIMENT,
                  group = interaction(EXPERIMENT,INDIVIDUAL)), alpha = 0.4) +
  scale_fill_manual(values = palette_exp, name = "Experiment")

# Facet by experiment
plot_nclones_sd_facets <- plot_nclones_sd +
  facet_grid(.~EXPERIMENT) +
  theme(strip.text = element_blank(), panel.spacing = unit(0.5, "cm"))

# Combine total and faceted plot
outplot_nclones <- gplot_grid_onelegend(plot_nclones_sd,plot_nclones_sd_facets,
                                        plot_height=plot_height, 
                                        plot_width=plot_width)

#------------------------------------------------------------------------------
# PLOT RAREFACTION CURVES OF P20
#------------------------------------------------------------------------------

# Plot mean lines only
plot_p_means <- ggplot(r %>% filter(METRIC == paste0("P", P))) +
  geom_line(aes(x=SAMPLE_SIZE, y=MEAN, colour = EXPERIMENT, 
                group = interaction(EXPERIMENT,INDIVIDUAL))) +
  scale_color_manual(values = palette_exp, name = "Experiment") +
  scale_x_continuous(breaks = seq(0, 10000, 2000), 
                     labels = function(x) x/1000, 
                     name = paste(x_title, "(×1000)")) +
  scale_y_continuous(breaks = seq(0,1,0.2),
                     labels = function(y) y * 100,
                     name = paste0("P", P, " (%)")) +
  theme_classic() + theme_base

# Add SD regions
plot_p_sd <- plot_p_means +
  geom_ribbon(aes(x=SAMPLE_SIZE, ymin=MEAN-SD, ymax = MEAN + SD, 
                  fill = EXPERIMENT,
                  group = interaction(EXPERIMENT,INDIVIDUAL)), alpha = 0.4) +
  scale_fill_manual(values = palette_exp, name = "Experiment")

# Facet by experiment
plot_p_sd_facets <- plot_p_sd +
  facet_grid(.~EXPERIMENT) +
  theme(strip.text = element_blank(), panel.spacing = unit(0.5, "cm"))

# Combine total and faceted plot
outplot_p <- gplot_grid_onelegend(plot_p_sd, plot_p_sd_facets,
                                  plot_height = plot_height, 
                                  plot_width = plot_width)


#------------------------------------------------------------------------------
# SAVE OUTPUT
#------------------------------------------------------------------------------

savefig(outplot_nclones, filename = paste0(outpath, "-counts"),
        height = plot_height, width = plot_width)
savefig(outplot_p, filename = paste0(outpath, "-p", P),
        height = plot_height, width = plot_width)