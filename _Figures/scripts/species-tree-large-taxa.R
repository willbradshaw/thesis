###############################################################################
## FIGURE                                                                    ##
## Small cladogram with Nfu/Xma/Ola/Gac IGH phenotypes                       ##
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

# Path to treefiles
tree_path <- "../_Data/trees/stree-cladogram.tre"

# Path to annotation data
annot_path <- "../_Data/species_data/species_tree_large.csv"


# Configure output
filename <- "species-tree-large-taxa"

#------------------------------------------------------------------------------
# Import and subset trees
#------------------------------------------------------------------------------

# Import trees
tree <- read.tree(tree_path)
tree$edge.length <- 1

# Import tree and annotations
annotation <- suppressMessages(read_csv(annot_path))

# Rearrange annotations to match tree (for clade labelling)
annotation <- arrange(annotation, match(tree_name, tree$tip.label))

#------------------------------------------------------------------------------
# Plot tree with taxon annotations
#------------------------------------------------------------------------------

# Basic tree topology and tip labels
treeplot <- suppressWarnings(revts(ggtree(tree)) + 
                               geom_rootedge(rootedge = 0.5))
t <- treeplot %<+% annotation
t <- t + geom_tiplab(aes(label = paste(Genus, Species), 
                         fontface = Highlight+3), size=6.5, offset=0.3,
                     family = font)

# Label Atherinomorpha
mrca_atherinomorpha <- ggtree::MRCA(tree, 
                                    which(annotation$Atherinomorpha == 1))
t <- t + geom_cladelabel(mrca_atherinomorpha, "Atherinomorpha",
                         offset = 18, offset.text = 0.75, fontsize=7,
                         extend = 0.2, barsize = 0.75, family = font,
                         angle = 270, hjust = 0.5)

# Label Cyprinidontiformes
mrca_cyprinidontiformes <- ggtree::MRCA(tree, 
                                      which(annotation$Cyprinodontiformes == 1))
t <- t + geom_cladelabel(mrca_cyprinidontiformes, "Cyprinidontiformes",
                         offset = 15, offset.text = 0.75, fontsize=7,
                         extend = 0.2, barsize = 0.75, family = font,
                         angle = 270, hjust = 0.5)

# Label Nothobranchiidae
mrca_nothobranchiidae <- ggtree::MRCA(tree, 
                                      which(annotation$Nothobranchiidae == 1))
t <- t + geom_cladelabel(mrca_nothobranchiidae, "Nothobranchiidae",
                         offset = 12, offset.text = 0.75, fontsize=7,
                         extend = 0.2, barsize = 0.75, family = font,
                         angle = 270, hjust = 0.5)

# Fix x-limit
t <- t + xlim(c(NA, 19))

#------------------------------------------------------------------------------
# Plot tree
#------------------------------------------------------------------------------

plot_height <- 30
plot_width <- 30

plt <- grid.grab()
savefig(plot = t, 
        filename = filename,
        height = plot_height, width = plot_width)