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
tree_path <- "../_Data/trees/stree-cladogram.tre"

# Path to annotation data
annot_path <- "../_Data/species_data/multispecies_cz_subclasses.csv"

# Configure output
filename <- "multispecies-cz-subclasses"

#------------------------------------------------------------------------------
# Import and subset tree topology
#------------------------------------------------------------------------------

# Import tree and annotations
annotation <- suppressMessages(read_csv(annot_path))
tree <- read.tree(tree_path)

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
tree$root.edge <- 0.5
tree$edge.length <- NA

# Basic tree topology
treeplot <- suppressWarnings(revts(ggtree(tree)) + geom_rootedge())
t <- treeplot

# Add line-breaks to species names
annotation$name <- paste0(substr(annotation$Genus, 1, 1), ". ", annotation$Species)

# Add annotations
label_size <- 7
t <- facet_plot(t, "Species", annotation, geom_text,
           aes(x=0, label=name), size = label_size, 
           fontface = "italic", family = font)
t <- facet_plot(t, "# IGHZ-A", annotation, geom_text,
                aes(x=0, label=IGHZA), size = label_size,
                family = font)
t <- facet_plot(t, "# IGHZ-B", annotation, geom_text,
                aes(x=0, label=IGHZB), size = label_size,
                family = font)
t <- facet_plot(t, "# IGHZ-C", annotation, geom_text,
                aes(x=0, label=IGHZC), size = label_size,
                family = font)

# Format facet strip labels

t <- t + theme(
  strip.text = element_text(size=fontsize_base * 1.6, face = "bold",
                            margin = margin(b=0.5, unit="cm"),
                            family = titlefont),
  strip.background = element_blank()
)

# Set viewport
plot_unit = "cm"
plot_height <- 20
plot_width <- plot_height * 1.2
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
# gt$widths[9] <- gt$widths[9] * 1.2
# gt$widths[17] <- gt$widths[17] * 1.2

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