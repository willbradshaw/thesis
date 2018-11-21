###############################################################################
## FIGURE                                                                    ##
## Nothobranchius furzeri VH families                                       ##
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

# Path to V sequences
vh_nt_path <- "../_Data/segments/nfu/nfu_vh_nt.fasta"

# Configure output
filename_comb <- "nfu-vh-families"
filename_tree <- "nfu-vh-families-tree"
filename_map <- "nfu-vh-families-map"

# Pairwise identity method
pid_type <- "PID2"

# Clustering method
clust_method <- "single"

# Clustering threshold for family membership
id_threshold <- 80

# Small font size for V labels
fontscale_small_label <- 0.5

#------------------------------------------------------------------------------
# PREPARE VH SEQUENCES AND PERFORM CLUSTERING
#------------------------------------------------------------------------------

# Import V sequences
vh_nt <- readDNAStringSet(vh_nt_path)
nv <- length(vh_nt)

# Remove terminal Ns
vh_nt <- DNAStringSet(sub("^N*", "", vh_nt))
vh_nt <- DNAStringSet(sub("N*$", "", vh_nt))

# Generate pairwise alignments and compute % ID
vh_id_nt <- sapply(seq(nv), function(subject)
  pid(pairwiseAlignment(vh_nt, vh_nt[subject]), 
      type = pid_type))
# Fix row and column names
rownames(vh_id_nt) <- names(vh_nt)
colnames(vh_id_nt) <- names(vh_nt)

# Convert PID into a distance measure
vh_dist_nt <- as.dist(100-vh_id_nt)

# Construct cluster dendrogram
vh_clust_nt <- hclust(vh_dist_nt, method = clust_method)

#------------------------------------------------------------------------------
# CONVERT TO PHYLOGENETIC TREE AND ADD FAMILY INFORMATION
#------------------------------------------------------------------------------

# Auxiliary functions (TODO: move to aux file?)
grepl_tips <- function(pattern, tree_object){
  tree_object$tip.label[grepl(pattern, tree_object$tip.label)]
}
get_mrca <- function(pattern, tree_object){
  nodes <- grepl_tips(pattern, tree_object)
  mrca <- ggtree::MRCA(tree_object, nodes)
  return(ifelse(is.null(mrca), NA, mrca))
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Convert dendrogram to phylogenetic tree
vh_phylo <- as.phylo(vh_clust_nt)
vh_tbl_tree_raw <- as_data_frame(vh_phylo)
vh_tbl_tree_raw[["locus"]] <- sub("IGH(\\d)V.*", "\\1", vh_tbl_tree_raw$label) %>%
  as.integer

# Group by family and count
vh_families <- vh_tbl_tree_raw %>% 
  filter(!is.na(label)) %>%
  mutate(family = sub("IGH\\dV(.*)-.*", "\\1", label)) %>%
  group_by(family) %>%
  count() %>%
  mutate(colour = colours[["GR"]]) %>%
  arrange(desc(n))

# Set colours for multi-segment families
n_large <- sum(vh_families$n > 1)
vh_families$colour[1:n_large] <- gg_color_hue(n_large)

# Add rows for non-family (grey) and non-annotated (white, transparent) entries
vh_families <- vh_families %>% ungroup() %>%
  add_row(family = "gr", n = 0, colour = colours[["GR"]]) %>%
  add_row(family = "wh", n = 0, colour = colours[["WH"]])

# Find and add root nodes for each family
vh_families <- vh_families %>%
  mutate(node = unlist(sapply(family, function(f) 
    get_mrca(paste0(f,"-"), vh_phylo))))

# Add family root data to tree tibble and propagate down
vh_tbl_tree_fam <- left_join(vh_tbl_tree_raw,
                             vh_families %>% select(-n, -colour),
                             by = "node")
for (n in seq(nrow(vh_families))){
  family <- vh_families$family[n]
  root <- vh_families$node[n]
  if (!is.na(root)){
    progeny <- offspring(vh_tbl_tree_fam, root) %>% pull(node)
    vh_tbl_tree_fam[vh_tbl_tree_fam$node %in% progeny,]$family <- family
  }
}

# Fill in blanks with solo family name (if available) or fallback
copy_fill <- !is.na(vh_tbl_tree_fam$label) & is.na(vh_tbl_tree_fam$family)
vh_tbl_tree_fam[copy_fill,]$family <- sub("IGH\\dV(.*)-.*", "\\1", 
                                          vh_tbl_tree_fam[copy_fill,]$label)
vh_tbl_tree_fam[is.na(vh_tbl_tree_fam$family),]$family <- "gr"

# Convert family to factor for palette ordering
vh_tbl_tree_fam$family <- factor(vh_tbl_tree_fam$family, 
                                 levels = vh_families$family)

# Correct branch length relative to original cluster object
vh_tbl_tree_fam$branch.length <- vh_tbl_tree_fam$branch.length / 
  (max(vh_tbl_tree_fam$branch.length, na.rm = TRUE) / max(vh_clust_nt$height))

#------------------------------------------------------------------------------
# PLOT VH DENDROGRAM (AS PHYLOGENETIC TREE)
#------------------------------------------------------------------------------

theme_tre <- theme_minimal() + theme_base +
  theme(legend.position = "top",
        legend.margin = margin(b=0,t=0.5, unit="cm"),
        plot.margin = margin(t=0, unit="cm"),
        axis.title.x = element_text(margin = margin(t=0, b=0.5, unit="cm")),
        axis.title.y = element_text(margin = margin(l = 2, t = 2, unit = "cm")))

tre <- revts(ggtree(as.treedata(vh_tbl_tree_fam), aes(colour=family), 
              ladderize=TRUE)) +
  scale_colour_manual(values = vh_families$colour, breaks = vh_families$family[1:n_large],
                      name = "VH family") +
  coord_flip() + 
  scale_x_reverse(labels = function(x) x + 100, name = "NT sequence identity (%)",
                  limits = c(5, NA)) +
  scale_y_continuous(breaks=NULL, name = "VH segment") +
  geom_tiplab(angle = 270, family = titlefont, size=2, 
              aes(fontface = as.integer(locus))) +
  geom_vline(xintercept = -20, colour = "red", linetype = "dashed") +
  guides(colour = guide_legend(nrow=1)) +
  theme_tre
tre

#------------------------------------------------------------------------------
# PREPARE DATA FOR HISTOGRAM
#------------------------------------------------------------------------------

# Determine order from tree data
vh_tree_order <- tre$data %>% arrange(y) %>% filter(!is.na(label)) %>% pull(label)

# Generate %identity table for VH segments
vh_tab <- vh_id_nt %>%
  melt(varnames = c("Query", "Subject"), value.name = "ID") %>%
    mutate(Close = ID >= id_threshold,
           Family = ifelse(sub("IGH\\dV(.*)-.*", "\\1", Query) == 
                             sub("IGH\\dV(.*)-.*", "\\1", Subject),
                           sub("IGH\\dV(.*)-.*", "\\1", Query), "wh")) %>%
  mutate(Query = factor(Query, levels = vh_tree_order),
         Subject = factor(Subject, levels = rev(vh_tree_order)),
         Family = factor(Family, levels = vh_families %>% filter(family != "gr") %>% pull(family)))

heatmap_label_colours <- levels(vh_tab$Query) %>% 
  sub("IGH\\dV(.*)-.*", "\\1", .) %>% 
  match(vh_families$family) %>%
  vh_families$colour[.]

heatmap_label_faces <- levels(vh_tab$Query) %>% 
  sub("IGH(\\d.*)V.*", "\\1", .) %>% 
  as.integer


#------------------------------------------------------------------------------
# PLOT VH FAMILY HEATMAP
#------------------------------------------------------------------------------

# Make identity heatmap
fontsize_hm_ticks <- fontsize_base/2
theme_hm_vh <- theme_classic() + theme_base +
  theme(axis.text.x = element_text(size=fontsize_hm_ticks, angle = 90),
        axis.text.y = element_text(size=fontsize_hm_ticks),
        legend.position = "top",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )

vh_heatmap <- ggplot(vh_tab) +
  geom_tile(aes(x=Query, y=Subject, fill=Family)) +
  geom_point(aes(x=Query, y=Subject, colour=Close), size=0.5) +
  scale_fill_manual(values = paste0(vh_families %>% filter(family != "gr") %>% 
                                      pull(colour), "66"),
                    breaks = vh_families$family[1:n_large],
                    name = "VH family") +
  scale_colour_manual(values = c("#FFFFFF00", "red"), guide=FALSE) +
  guides(fill = guide_legend(nrow=1)) +
  theme_hm_vh + coord_fixed() +
  theme(axis.text.x = element_text(colour = heatmap_label_colours,
                                   face = heatmap_label_faces),
        axis.text.y = element_text(colour = rev(heatmap_label_colours),
                                   face = rev(heatmap_label_faces)))

#------------------------------------------------------------------------------
# STACK PLOTS
#------------------------------------------------------------------------------

rel_height_tree <- 0.6
plt <- plot_grid(tre, vh_heatmap, ncol = 1, labels="AUTO", 
                 label_fontfamily = titlefont, label_fontface = "plain",
                 label_size = fontsize_base * fontscale_label,
                 rel_heights = c(rel_height_tree, 1))



#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plot_height = 20
plot_ratio_tree = 1
plot_ratio_map = 1/1.05
plot_ratio_comb = 1/(1.05 + rel_height_tree)

savefig(plot = tre, filename = filename_tree, 
        height = plot_height, ratio = plot_ratio_tree)
savefig(plot = vh_heatmap, filename = filename_map,
        height = plot_height, ratio = plot_ratio_map)
savefig(plot = plt, filename = filename_comb,
        height = plot_height*1.5, ratio = plot_ratio_comb)