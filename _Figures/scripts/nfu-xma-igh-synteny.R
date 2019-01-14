###############################################################################
## FIGURE                                                                    ##
## Synteny plot of chromosomes containing IGH locu in Xma and Nfu            ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")
source("aux/ggplot2.R")
source("aux/trees.R")

# Path to links files
nfu_xma_links_path <- "../_Data/synteny/links_xmac_blocks.txt"

# Specify chromosomes to match
species <- c("nfu", "xma")
n_nfu <- which(species == "nfu")
n_xma <- which(species == "xma")
chrs <- tibble(
  species = species,
  name = c("chr6", "group16"),
  length = c(61955777, 25094516),
  igh_start = c(43205521 + 98762 - 1, 9384676 + 98855 - 1),
  igh_end = c(43205521 + 253166 - 1, 9384676 + 361512 - 1)
)

# Output path
filename <- "nfu-xma-igh-synteny"

#------------------------------------------------------------------------------
# PREPARE SYNTENY PLOT BETWEEN CHROMOSOMES
#------------------------------------------------------------------------------

# Load links table
colnames_links <- paste(rep(species, each=3), c("chr", "start", "end"), sep="_")
links <- suppressMessages(read_delim(nfu_xma_links_path, delim = " ", 
                    col_names = c(colnames_links, "comment"))) %>%
  filter(nfu_chr == chrs$name[n_nfu], xma_chr == chrs$name[n_xma]) %>%
  mutate(colour = colours[["LG"]])
  
# Define DNA segments to plot (whole chromosome + IGH locus)
dna_segs_df <- lapply(species, function(s) chrs %>% filter(species == s) %>%
                        transmute(name=name, start=1, end=length, strand=1, 
                                  col = colours[[s]], fill = colours[[s]],
                                  gene_type = "lines"))
igh_segs_df <- lapply(species, function(s) chrs %>% filter(species == s) %>%
                        transmute(name="IGH", start=igh_start, end=igh_end, 
                                  strand=1, gene_type = "blocks",
                                  col = colours[[s]], fill = colours[[s]]))

dna_segs <- lapply(seq(length(dna_segs_df)), function(n)
  dna_seg(bind_rows(dna_segs_df[[n]], igh_segs_df[[n]])))

# Define link comparisons to plot
comparisons_df <- links %>%
  transmute(start1 = nfu_start, end1 = nfu_end, start2 = xma_start, 
            end2 = xma_end, col=colour)
comparisons <- comparisons_df %>% as.comparison() %>% list

# Define ranges to plot
xlims <- lapply(dna_segs_df, function(x) c(x$start, x$end))

# Add names
seg_names <- c(expression(italic(N.~furzeri)~~~chromosome~6),
               expression(italic(X.~maculatus)~~~chromosome~16))
names(dna_segs) <- seg_names
names(xlims) <- seg_names

# Create IGH annotation objects
annotations_df <- lapply(species, function(s) chrs %>% filter(species == s) %>%
                           transmute(x1 = igh_start, x2 = igh_end, text = "IGH", 
                                     color = colours[[s]]))
annotations <- lapply(annotations_df, function(a) as.annotation(a))

#------------------------------------------------------------------------------
# CREATE AND REFINE SYNTENY PLOT
#------------------------------------------------------------------------------

# Define basic parameters
plot_height <- 5
plot_width <- 25
plot_unit <- "cm"

# Create viewport
layout <- grid.layout(
  ncol = 1, nrow = 1,
  heights = unit(plot_height, plot_unit),
  widths = unit(plot_width, plot_unit)
)
v <- viewport(layout = layout, name="main")
grid.newpage()
pushViewport(v)

# Create plot
plot_gene_map(
  dna_segs = dna_segs,
  comparisons = comparisons,
  tree = NULL, tree_width = NULL, tree_branch_labels_cex = NULL, tree_scale = FALSE, legend = NULL,
  annotations = NULL,
  annotations_height = 3,
  annotations_cex = 1.5,
  seg_plots = NULL, seg_plot_height = 1, seg_plot_height_unit = "lines", seg_plot_yaxis = 3, seg_plot_yaxis_cex = scale_cex,
  xlims = xlims,
  offsets = NULL,
  minimum_gap_size = 0.05, fixed_gap_length = FALSE, 
  limit_to_longest_dna_seg = FALSE,
  main = NULL, main_pos = "centre", dna_seg_labels = NULL,
  dna_seg_label_cex = 1, dna_seg_label_col = "black",
  gene_type = NULL,
  arrow_head_len = 200, 
  dna_seg_line = FALSE,
  scale = TRUE, dna_seg_scale = FALSE, n_scale_ticks = 7, scale_cex = 0.3,
  global_color_scheme = c("auto", "auto", "blue_red", 0.5), override_color_schemes = FALSE,
  plot_new = FALSE,
  debug = 0
)
upViewport(0)

# Edit appearance
downViewport("main")
for (i in 1:length(seg_names)){
  grid.edit(paste("label", i, sep="."), label=seg_names[i], gp = gpar(fontfamily = font, fontsize=fontsize_base)) # ...
  grid.edit(paste("seg", i, chrs$name[i], sep="."), gp = gpar(lwd=5))
}
# ...
upViewport(0)

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plt <- grid.grab()
savefig(plot = plt, filename = filename, dpi = 640,
        height = plot_height, width = plot_width)