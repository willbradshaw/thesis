###############################################################################
## FIGURE                                                                    ##
## Phylogenetic tree of multispecies CZ sequences                            ##
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
source("aux/trees.R")

# Path to treefiles
tre_nt_path <- "../_Data/trees/multispecies-cz-nt.tre"

# Configure paremeters
min_support <- 65 # All lower supports considered unsupported
species_codes <- c("xma", "pre", "pfo", "fhe", "cto", "kma",
                   "cva", "ppl", "nfu")
species_names <- c("X. maculatus", "P. reticulata", "P. formosa",
                   "F. heteroclitus", "C. toddi", "K. marmoratus",
                   "C. variegatus", "P. playfairii", "N. furzeri")
rootedge_length <- 0.1
treeline_width <- 1.1
tiplab_offset <- 0.035
zclade_patterns <- c("P. formosa IGHZ1|K. marmoratus IGHZ1",
                     "P. formosa IGHZ2|C. toddi IGHZ4",
                     "P. formosa IGHZ3|K. marmoratus IGHZ2")
zclade_labels <- LETTERS[1:length(zclade_patterns)]
zclade_colours <- gg_color_hue(3)

# Configure output
outfile_nt <- "multispecies-cz-tree"

#------------------------------------------------------------------------------
# IMPORT TREE AND CORRECT TIP LABELS
#------------------------------------------------------------------------------

tre_nt <- read.tree(tre_nt_path)

# Replace species codes with names
for (n in 1:length(species_codes)){
  tre_nt$tip.label <- gsub(paste0(species_codes[n], "_"),
                           paste0(species_names[n], " "),
                           tre_nt$tip.label)
}

# Remove extra annotations
tre_nt$tip.label <- gsub("_.*", "", tre_nt$tip.label)

#------------------------------------------------------------------------------
# COLLAPSE UNSUPPORTED NODES
#------------------------------------------------------------------------------

tre_nt <- as.polytomy(tre_nt,'node.label',
                   function(x) as.numeric(x) < min_support)

#------------------------------------------------------------------------------
# ADD ISOTYPE AND SUBCLADE ANNOTATION
#------------------------------------------------------------------------------

# Get MRCAs of IGHZ and IGHM (outgroup)
mrca_z <- get_mrca("IGHZ", tre_nt)
mrca_m <- get_mrca("IGHM", tre_nt)

# Set isotype status
nodes_z <- c(mrca_z, offspring(tre_nt, mrca_z))
nodes_m <- c(mrca_m, offspring(tre_nt, mrca_m))

# Convert into tree table and add isotype status
tab <- as_tibble(tre_nt) %>% mutate(isotype = "GR", zclade = NA)
tab$isotype[tab$node %in% nodes_z] <- "CZ"
tab$isotype[tab$node %in% nodes_m] <- "CM"
class(tab) <- c("tbl_tree", class(tab))

# Determine major IGHZ clades
zclade_mrcas <- sapply(zclade_patterns, function(p) get_mrca(p, tre_nt))
zclade_nodes <- lapply(zclade_mrcas, 
                       function(m) c(m, offspring(tab, m) %>% pull(node)))
# Exclude zclades with exactly one taxon
zclade_keep <- sapply(zclade_nodes, length) > 1
zclade_mrcas <- zclade_mrcas[zclade_keep]
zclade_nodes <- zclade_nodes[zclade_keep]

# Apply zclade labels to table
for (n in 1:length(zclade_nodes)){
  tab$zclade[tab$node %in% zclade_nodes[[n]]] <- n 
}

#------------------------------------------------------------------------------
# PLOT AND ANNOTATE TREE
#------------------------------------------------------------------------------
treedata_nt <- as.treedata(tab)
treedata_nt@phylo$root.edge <- rootedge_length

# Make basic plot and colour by isotype
treeplot <- ggtree(treedata_nt, aes(colour = isotype), size = treeline_width) +
  geom_rootedge(colour = colours[["GR"]], size = treeline_width) +
  geom_tiplab(size = 6.5, offset = tiplab_offset, family = font) +
  geom_nodelab(nudge_x = 0.003, hjust = 0, size = 3) + 
  scale_colour_manual(values = unlist(colours[tab$isotype]))

# Add zclade labels
zclade_offsets <- rep(0.5,3)
for (n in 1:length(zclade_labels)){
  treeplot <- treeplot + 
    geom_cladelabel(zclade_mrcas[n], zclade_labels[n], family = font,
                    offset = zclade_offsets[n], offset.text = 0.05,
                    fontsize = 8, barsize = 0.75, hjust = 0.5,
                    extend = 0.3, align = TRUE, color = zclade_colours[n])
}
                                         
#------------------------------------------------------------------------------
# SAVE TREE PLOT
#------------------------------------------------------------------------------
plot_height = 25
plot_ratio = 1

savefig(treeplot, outfile_nt, height = plot_height, ratio = plot_ratio)