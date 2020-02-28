###############################################################################
## FIGURE                                                                    ##
## Species cladogram                                                         ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# Import and subset trees
#------------------------------------------------------------------------------

# Import trees
tree <- read.tree(inpath_tree)
tree$edge.length <- 1

# Import tree and annotations
annotation <- suppressMessages(read_csv(inpath_annot))

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
mrca_atherinomorpha <- MRCA(tree, 
                                    which(annotation$Atherinomorpha == 1))
t <- t + geom_cladelabel(mrca_atherinomorpha, "Atherinomorpha",
                         offset = 18, offset.text = 0.75, fontsize=7,
                         extend = 0.2, barsize = 0.75, family = font,
                         angle = 270, hjust = 0.5)

# Label Cyprinidontiformes
mrca_cyprinidontiformes <- MRCA(tree, 
                                      which(annotation$Cyprinodontiformes == 1))
t <- t + geom_cladelabel(mrca_cyprinidontiformes, "Cyprinidontiformes",
                         offset = 15, offset.text = 0.75, fontsize=7,
                         extend = 0.2, barsize = 0.75, family = font,
                         angle = 270, hjust = 0.5)

# Label Nothobranchiidae
mrca_nothobranchiidae <- MRCA(tree, 
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
ggsave(plot = t, filename = outpath, device = "svg", units = "cm",
       height=plot_height, width = plot_width)
