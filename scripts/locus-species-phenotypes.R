###############################################################################
## FIGURE                                                                    ##
## Small cladogram with Nfu/Xma/Ola/Gac IGH phenotypes                       ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# Import and subset tree topology
#------------------------------------------------------------------------------

write_log("Importing data...", newline = FALSE)
# Import tree and annotations
annotation <- suppressMessages(read_csv(inpath_annot))
tree <- read.nexus(inpath_tree, tree.names = "Teleostei")
log_done(newline=TRUE)

# Prune tree to desired species
write_log("Pruning tree to desired species...", newline = FALSE)
prune_tree <- function(tree, nodes_keep){
  drop.tip(tree, tree$tip.label[-nodes_keep])
}
nodes_keep <- match(annotation$tree_name, tree$tip.label)
nodes_keep <- nodes_keep[!is.na(nodes_keep)]
tree <- prune_tree(tree, nodes_keep)
log_done(newline=TRUE)

#------------------------------------------------------------------------------
# Plot tree with annotation table
#------------------------------------------------------------------------------

# Basic tree topology
write_log("Plotting tree...", newline = FALSE)
treeplot <- revts(ggtree(tree)) + geom_rootedge(rootedge=0.2)
t <- treeplot
log_done(newline=TRUE)

write_log("Adding annotations...")
# Add line-breaks to species names
annotation$common_name <- gsub(" ", "\n", annotation$common_name)
annotation$scientific_name <- gsub(" ", "\n", annotation$scientific_name)

# Add annotations
label_size <- 7
t <- facet_plot(t, "Species", annotation, geom_text,
           aes(x=0, label=scientific_name), size = label_size, 
           fontface = "italic", family = font)
write_log("    1")
t <- facet_plot(t, "Common\nname", annotation, geom_text,
                aes(x=0, label=common_name), size = label_size,
                family = font)
write_log("    2")
t <- facet_plot(t, "IGHZ?", annotation, geom_text,
                aes(x=0, label=c("No", "Yes")[has_ighz+1], 
                    fontface = 2-has_ighz), size = label_size,
                family = font)
write_log("    3")
t <- facet_plot(t, "Inverted\nsublocus?", annotation, geom_text,
                aes(x=0, label=c("No", "Yes")[has_inverted_sublocus+1],
                    fontface = has_inverted_sublocus+1), size = label_size,
                family = font)
write_log("    4")
t <- facet_plot(t, "# IGHM-TM\nexons", annotation, geom_text,
                aes(x=0, label=n_ighm_tm_exons, fontface=6-n_ighm_tm_exons),
                size = label_size, family = font)
write_log("    5")
t <- facet_plot(t, "# CD(2,3,4)\nduplications", annotation, geom_text,
                aes(x=0, label=cd234_duplications, 
                    fontface = cd234_duplications), size = label_size,
                family = font)
write_log("    6")
log_done(newline=TRUE)

# Format facet strip labels

t <- t + theme(
  strip.text = element_text(size=fontsize_base * 1.6, face = "bold",
                            margin = margin(b=1.2, unit="cm"),
                            family = titlefont),
  strip.background = element_blank()
)

write_log("Writing to grid...", newline = FALSE)
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
log_done(newline=TRUE)

#------------------------------------------------------------------------------
# Plot tree with annotation table
#------------------------------------------------------------------------------

write_log("Saving output...", newline = FALSE)
plt <- grid.grab()
ggsave(plot = plt, filename = outpath, device = "svg", units = "cm",
       height=plot_height, width = plot_width)
log_done(newline=TRUE)
