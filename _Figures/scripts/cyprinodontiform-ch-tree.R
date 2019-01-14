###############################################################################
## FIGURE                                                                    ##
## Phylogenetic tree of cyprinodontiform CH domains                          ##
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
tre_nt_path <- "../_Data/trees/ch-nt.tre"

# Configure parameters
exon_patterns <- c("IGH.?M.?-1", "IGH.?M.?-2", "IGH.?M.?-3", "IGH.?M.?-4",
                 "IGH.?Z.?-.*1_?$", "IGH.?Z.?-.*2_?$", "IGH.?Z.?-.*3_?$", "IGH.?Z.?-.*4_?$",
                 "IGH.?D.?-1", "IGH.?D.?-2", "IGH.?D.?-3", "IGH.?D.?-4",
                 "IGH.?D.?-5", "IGH.?D.?-6", "IGH.?D.?-7")
exon_labels <- c("CM1", "CM2", "CM3", "CM4",
                 "CZ1", "CZ2", "CZ3", "CZ4",
                 "CD1", "CD2", "CD3", "CD4",
                 "CD5", "CD6", "CD7")
exon_fills <- sapply(exon_labels, function(x) colours[[gsub("\\d", "", x)]])
exon_offsets <- c(-0.04, 0.01, 0.02, 0.055,
                  -0.01, -.05, 0, -0.01,
                  0.04, -0.04, -0.03, 0.02,
                  0.0556, 0.03, 0.04)

# Configure output
outfile_nt <- "ch-tree-all"
exon_tree_path_prefix <- "../_Data/trees/ch-tree-exon"
outfile_exon_prefix <- "ch-tree-exon"
outfile_exon <- paste(outfile_exon_prefix, tolower(exon_labels), sep="-")
names(outfile_exon) <- exon_labels

# Small font size for V labels
fontscale_small_label <- 0.5

#------------------------------------------------------------------------------
# IMPORT ALL-EXON TREE
#------------------------------------------------------------------------------

tre_nt <- read.tree(tre_nt_path)

#------------------------------------------------------------------------------
# DETERMINE GROUP (EXON) MRCAS
#------------------------------------------------------------------------------

exon_mrcas_nt <- sapply(exon_patterns, function(x) get_mrca(x, tre_nt))
names(exon_mrcas_nt) <- exon_labels

#------------------------------------------------------------------------------
# MAKE ALL-EXON TREE
#------------------------------------------------------------------------------

tree_layout <- "daylight"
fontscale_treelabel <- 0.5

ch_parent_tree <- function(tree, layout, mrcas, fills,
                           marker_size=4, marker_shape=23){
  # Make base tree
  g_base <- ggtree(tree, layout=layout)
  g <- g_base
  # Collapse nodes
  for (n in seq(length(mrcas))){
    g <- g + 
      geom_cladelabel2(node=mrcas[n], label=names(mrcas)[n], family = font,
                       barsize = 0, offset = exon_offsets[n], hjust = 0.5,
                       fontsize = fontsize_base * fontscale_treelabel) + 
      geom_hilight_encircle(node=mrcas[n], fill = fills[n], expand = 0.01,
                            s_shape = 1)
  }
  return(g)
}

allplot_nt <- ch_parent_tree(tre_nt, tree_layout, exon_mrcas_nt, exon_fills)

#------------------------------------------------------------------------------
# SAVE ALL-EXON TREE
#------------------------------------------------------------------------------
# 
plot_height = 25
plot_ratio = 1

savefig(plot = allplot_nt, filename = outfile_nt,
        height = plot_height, ratio = plot_ratio)

# #------------------------------------------------------------------------------
# # CUT TREE TO OBTAIN EXON TREES
# #------------------------------------------------------------------------------
# 
# strees <- subtrees(tre_nt)
# exon_trees <- list()
# 
# for (n in seq(length(exon_patterns))){
#   # Get list of matching exons
#   exons <- sort(tre_nt$tip.label[grepl(exon_patterns[n], tre_nt$tip.label)])
#   # Get subtree with matching tip list
#   all_in <- sapply(strees, function(tree) 
#     all(exons %in% tree$tip.label) && all(tree$tip.label %in% exons))
#   # Assign to exon_trees
#   exon_trees[[exon_labels[n]]] <- strees[all_in][[1]]
# }
# 
# #------------------------------------------------------------------------------
# # SAVE EXON TREE DATA
# #------------------------------------------------------------------------------
# 
# # for (t in names(exon_trees)){
# #   tpath <- paste(exon_tree_path_prefix, tolower(t), sep = "-")
# #   tpath <- paste0(tpath, ".tre")
# #   write.tree(exon_trees[[t]], tpath)
# # }
# 
# #------------------------------------------------------------------------------
# # VISUALISE EXON TREES
# #------------------------------------------------------------------------------
# 
# # Extract support values into new variable
# min_support <- 50 # All lower supports considered unsupported
# exon_trees_support <- list()
# for (n in seq(length(exon_patterns))){
#   lab <- exon_labels[n]
#   tab <- as_tibble(exon_trees[[lab]])
#   tab$support <- as.numeric(tab$label)
#   tab$support <- ifelse(is.na(tab$support), 100, tab$support)
#   tab$support <- pmax(min_support, tab$support)
#   exon_trees_support[[lab]] <- as.treedata(tab)
# }
# 
# # Configure tree objects
# branch_cols <- list(
#   "CM" = colours[["CM"]],
#   "CD" = MixColor(colours[["CD"]], "white", 0.7),
#   "CZ" = colours[["CZ"]]
# )
# exon_xlims <- c(0.55, 0.95, 0.75, 0.5,
#                 0.65, 0.7, 0.82, 0.6, 
#                 0.8, 0.9, 0.55, 0.35, 0.45, 0.95, 0.65)
# names(exon_xlims) <- exon_labels
# exon_trees_vis <- list()
# col_low <- "black"
# branch_size <- 1
# tiplab_size <- 4
# tiplab_offset <- 0.005
# rootedge_length <- 0.05
# 
# # Generate tree objects
# for (n in names(exon_trees_support)){
#   exon_type <- sub("(C.).", "\\1", n)
#   col_high <- branch_cols[[exon_type]]
#   root_support <- as_tibble(tre_nt) %>% filter(node == exon_mrcas_nt[n]) %>% 
#     pull(label) %>% as.numeric()
#   root_colour <- MixColor(col_low, col_high, 1-root_support/100)
#   t <- ggtree(exon_trees_support[[n]], aes(colour = support), size=branch_size) +
#     geom_tiplab(family = titlefont, size = tiplab_size, offset = tiplab_offset) + 
#     geom_rootedge(rootedge = rootedge_length, colour = root_colour, size = branch_size) + 
#     scale_color_continuous(low=col_low, high=col_high, limits=c(min_support,100)) +
#     xlim_tree(exon_xlims[n]) +
#     theme_base + theme_tre + theme(legend.position = "bottom", 
#                                    legend.key.width = unit(1, units="cm"))
#   exon_trees_vis[[n]] <- t
# }
# 
# #------------------------------------------------------------------------------
# # SAVE EXON TREES
# #------------------------------------------------------------------------------
# exon_plot_height = 25
# exon_plot_ratio = 1
# 
# for (l in exon_labels){
#   savefig(plot = exon_trees_vis[[l]],
#           filename = outfile_exon[l],
#           height = exon_plot_height, ratio = exon_plot_ratio)
# }
# 
# # TODO: Make plot grid of multiple exon trees (e.g. all CM, all CZ)