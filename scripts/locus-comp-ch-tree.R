###############################################################################
## FIGURE                                                                    ##
## Phylogenetic tree of cyprinodontiform CH domains                          ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------
# IMPORT ALL-EXON TREE
#------------------------------------------------------------------------------

tre_nt <- read.tree(inpath)

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

ggsave(plot = allplot_nt, filename = outpath, device = "svg", units = "cm",
       height=plot_height, width=plot_height*plot_ratio)
