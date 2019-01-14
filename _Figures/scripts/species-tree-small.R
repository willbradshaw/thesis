###############################################################################
## FIGURE                                                                    ##
## Small cladogram with Nfu/Xma/Ola/Gac IGH phenotypes                       ##
###############################################################################
# TODO: Add dre?

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

# Path to treefile
tree_path <- "../_Data/trees/Hughes18_31Calibrations_Chronogram_nexus.tre"

# Path to annotation data
annot_path <- "../_Data/species_data/igh_tree_small.csv"

# Paths to pictures
# ...

# Configure output
filename <- "species-tree-small"

#------------------------------------------------------------------------------
# Import and subset tree topology
#------------------------------------------------------------------------------
# TODO: Include takifugu?

# Import tree and annotations
annotation <- suppressMessages(read_csv(annot_path))
tree <- read.nexus(tree_path, tree.names = "Teleostei")

# Prune tree to desired species
prune_tree <- function(tree, nodes_keep){
  drop.tip(tree, tree$tip.label[-nodes_keep])
}
nodes_keep <- match(annotation$tree_name, tree$tip.label)
nodes_keep <- nodes_keep[!is.na(nodes_keep)]
tree <- prune_tree(tree, nodes_keep)

#------------------------------------------------------------------------------
# Plot tree with annotation table
#------------------------------------------------------------------------------

# Basic tree topology
treeplot <- revts(ggtree(tree)) + geom_rootedge(rootedge=0.2)
t <- treeplot

# Add line-breaks to species names
annotation$common_name <- gsub(" ", "\n", annotation$common_name)
annotation$scientific_name <- gsub(" ", "\n", annotation$scientific_name)

# Add annotations
label_size <- 7
t <- facet_plot(t, "Species", annotation, geom_text,
           aes(x=0, label=scientific_name), size = label_size, 
           fontface = "italic", family = font)
t <- facet_plot(t, "Common\nname", annotation, geom_text,
                aes(x=0, label=common_name), size = label_size,
                family = font)
t <- facet_plot(t, "IGHZ?", annotation, geom_text,
                aes(x=0, label=c("No", "Yes")[has_ighz+1], 
                    fontface = 2-has_ighz), size = label_size,
                family = font)
t <- facet_plot(t, "Inverted\nsublocus?", annotation, geom_text,
                aes(x=0, label=c("No", "Yes")[has_inverted_sublocus+1],
                    fontface = has_inverted_sublocus+1), size = label_size,
                family = font)
t <- facet_plot(t, "# IGHM-TM\nexons", annotation, geom_text,
                aes(x=0, label=n_ighm_tm_exons, fontface=6-n_ighm_tm_exons),
                size = label_size, family = font)
t <- facet_plot(t, "# CD(2,3,4)\nduplications", annotation, geom_text,
                aes(x=0, label=cd234_duplications, 
                    fontface = cd234_duplications), size = label_size,
                family = font)

# Format facet strip labels

t <- t + theme(
  strip.text = element_text(size=fontsize_base * 1.6, face = "bold",
                            margin = margin(b=1.2, unit="cm"),
                            family = titlefont),
  strip.background = element_blank()
)

# Set viewport
plot_unit = "cm"
plot_height <- 20
plot_width <- plot_height * 2
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

# Adjust facet widths
gt = ggplot_gtable(ggplot_build(t))
gt$widths[5] <- gt$widths[5] * 2
gt$widths[7] <- gt$widths[7] * 1.5
gt$widths[9] <- gt$widths[9] * 1.2
gt$widths[17] <- gt$widths[17] * 1.2

# Draw to viewport
grid.draw(gt)
popViewport(1)

#------------------------------------------------------------------------------
# Plot tree with annotation table
#------------------------------------------------------------------------------

plt <- grid.grab()
savefig(plot = plt, 
        filename = filename,
        height = plot_height, width = plot_width)