###############################################################################
## FIGURE                                                                    ##
## Phylogenetic tree of N. furzeri and X. maculatus VH families              ##
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
outgroup_pattern <- "TRB"

#------------------------------------------------------------------------------
# IMPORT PRE-ROOTED PHYLOGENETIC TREES
#------------------------------------------------------------------------------

tre_nt <- read.tree(inpath)

#------------------------------------------------------------------------------
# DETERMINE SPECIES AND FAMILY INFORMATION
#------------------------------------------------------------------------------

annotate_vh_tree <- function(tre, pattern){
  # Add species and family annotations to a VH tree
  tre_tab <- as_tibble(tre) %>% label_tips(TRUE) %>%
    column_from_label("species", "_.*", "") %>% 
    column_from_label("family", "(.*)_.*IGH.?V(.*)-.*", "\\1_\\2")
  tre_tab[grepl(outgroup_pattern, tre_tab[["label"]]),][["family"]] <- pattern
  tre_tab <- tre_tab %>% propagate_column("species") %>% 
    propagate_column("family")
  return(as.treedata(tre_tab))
}

treedata_nt <- annotate_vh_tree(tre_nt, outgroup_pattern)

#------------------------------------------------------------------------------
# IDENTIFY PARENT NODES FOR HOMOLOGOUS FAMILIES
#------------------------------------------------------------------------------
pattern_group1 <- "xma_IGHV02|nfu_IGH.V1"
pattern_group2 <- "xma_IGHV03|nfu_IGH.V[24]"

group1_mrca_nt <- get_mrca(pattern_group1, treedata_nt@phylo)
group2_mrca_nt <- get_mrca(pattern_group2, treedata_nt@phylo)
# Identify outgroup node
outgroup_nt <- check_monophyly(outgroup_pattern, treedata_nt@phylo)

# TODO: Add and refine group highlighting on tree plot

#------------------------------------------------------------------------------
# PLOT VH DENDROGRAM (AS PHYLOGENETIC TREE)
#------------------------------------------------------------------------------

n_spc <- 4

theme_tre <- theme_tree2() + theme_base +
  theme(legend.position = "top",
        legend.margin = margin(b=0,t=0.5, unit="cm"),
        plot.margin = margin(t=0, unit="cm"),
        axis.title.x = element_text(margin = margin(t=0, b=0.5, unit="cm")),
        axis.title.y = element_text(margin = margin(l = 2, t = 2, unit = "cm")))

treeplot_vh <- function(treedata, outgroup_node, hilight_nodes, hilight_fills,
                        hilight_extendto=2, hilight_alpha=0.1){
  # Make tree
  g <- ggtree(treedata, aes(colour=species), ladderize = TRUE) %>%
    collapse(node=outgroup_node) +
    geom_point2(aes(subset=(node == outgroup_node)), size=4, shape=23, 
                fill=gg_color_hue(n_spc)[1]) +
    geom_rootedge(0.2, colour="#888888") +
    xlim_tree(hilight_extendto)
    theme_tre
  # Add annotations
  for (n in seq(length(hilight_nodes))){
    g <- g + geom_hilight(node=hilight_nodes[n], extendto=hilight_extendto, 
                          alpha=hilight_alpha, fill=hilight_fills[n])
  }
  g <- g + geom_tiplab(family = titlefont, size = 1.5)
  # Put root at top
  root <- get_mrca(".*", treedata@phylo)
  top_children <- child(treedata@phylo, root)
  return(flip(g, top_children[1], top_children[2]))
}

treeplot_nt <- treeplot_vh(treedata_nt, outgroup_nt,
                           c(group1_mrca_nt, group2_mrca_nt),
                           c(colours[["HL1"]], colours[["HL2"]]), 2)

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plot_height = 20
plot_ratio = 0.67

ggsave(plot = treeplot_nt, filename = outpath, device = "svg", units = "cm",
       height=plot_height, width = plot_height*plot_ratio)
