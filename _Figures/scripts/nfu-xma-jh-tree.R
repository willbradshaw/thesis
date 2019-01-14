###############################################################################
## FIGURE                                                                    ##
## Phylogenetic tree of N. furzeri and X. maculatus VH families              ##
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
tre_nt_path <- "../_Data/trees/nfu-xma-jh-nt.tre"
tre_aa_path <- "../_Data/trees/nfu-xma-jh-aa.tre"

# Configure parameters
outgroup_pattern <- "TRB"

# Configure output
outfile_nt <- "nfu-xma-jh-tree-nt"
outfile_aa <- "nfu-xma-jh-tree-aa"

# Small font size for J labels
fontscale_small_label <- 0.5

internal_numeric <- TRUE

#------------------------------------------------------------------------------
# IMPORT PRE-ROOTED PHYLOGENETIC TREES
#------------------------------------------------------------------------------

tre_nt <- read.tree(tre_nt_path)
tre_aa <- read.tree(tre_aa_path)

#------------------------------------------------------------------------------
# DETERMINE SPECIES AND FAMILY INFORMATION
#------------------------------------------------------------------------------

annotate_jh_tree <- function(tre, pattern){
  # Add species and family annotations to a VH tree
  as_tibble(tre) %>% label_tips(internal_numeric) %>%
    column_from_label("species", "_.*", "") %>%
    propagate_column("species") %>% as.treedata
}

treedata_nt <- annotate_jh_tree(tre_nt, outgroup_pattern)
treedata_aa <- annotate_jh_tree(tre_aa, outgroup_pattern)

#------------------------------------------------------------------------------
# IDENTIFY SPECIAL NODES FOR ANNOTATION
#------------------------------------------------------------------------------

# Identify outgroup node
outgroup_nt <- check_monophyly(outgroup_pattern, treedata_nt@phylo)
outgroup_aa <- check_monophyly(outgroup_pattern, treedata_aa@phylo)

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

treeplot_jh <- function(treedata, outgroup_node){
  # Make tree
  g <- ggtree(treedata, aes(colour=species), ladderize = TRUE) %>%
    collapse(node=outgroup_node) +
    geom_point2(aes(subset=(node == outgroup_node)), size=4, shape=23, 
                fill=gg_color_hue(n_spc)[1]) +
    geom_rootedge(0.2, colour="#888888") +
    theme_tre
  # Add annotations
  g <- g + geom_tiplab(family = titlefont, size = 1.5)
  # Put root at top
  root <- get_mrca(".*", treedata@phylo)
  top_children <- child(treedata@phylo, root)
  return(flip(g, top_children[1], top_children[2]))
}

treeplot_nt <- treeplot_jh(treedata_nt, outgroup_nt)
treeplot_aa <- treeplot_jh(treedata_aa, outgroup_aa)

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plot_height = 15
plot_ratio = 1

savefig(plot = treeplot_nt, filename = outfile_nt,
        height = plot_height, ratio = plot_ratio)
savefig(plot = treeplot_aa, filename = outfile_aa,
        height = plot_height, ratio = plot_ratio)