###############################################################################
## FIGURE                                                                    ##
## Large cladogram with IGHZ annotations                                     ##
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

# Import tree and annotations
annotation <- suppressMessages(read_csv(inpath_annot))

# Rearrange annotations to match tree (for clade labelling)
annotation <- arrange(annotation, match(tree_name, tree$tip.label))

#------------------------------------------------------------------------------
# Infer IGHZ state for internal nodes
#------------------------------------------------------------------------------

# Convert tree to tibble and add annotations
tree_tab <- full_join(as_tibble(tree),
                      annotation %>% rename(label = tree_name),
                      by = "label")

# Propagate IGHZ status up from tips
ighz_status_unknown <- tree_tab %>% filter(is.na(IGHZ)) %>% pull(node)
while (length(ighz_status_unknown > 0)){
  # Find and count offspring for each node
  progeny <- lapply(ighz_status_unknown, function(n) offspring(tree_tab, n))
  n_progeny <- sapply(progeny, nrow)
  # Pick a parent with (joint) lowest number of offspring
  m <- sample(which(n_progeny == min(n_progeny)), 1)
  # Assign IGHZ status
  progeny_ighz <- progeny[[m]]$IGHZ
  if (any(is.na(progeny_ighz))) { ighz_status <- NA
  } else if (all(progeny_ighz == progeny_ighz[1])){ ighz_status <- progeny_ighz[1]
  } else { ighz_status <- "MIXED" }
  tree_tab[tree_tab$node == ighz_status_unknown[m],]$IGHZ <- ighz_status
  # Recompute unknown node list
  ighz_status_unknown <- tree_tab %>% filter(is.na(IGHZ)) %>% pull(node)
}

# Resolve mixed nodes by propagating down from root
ighz_status_root <- 1
root_id <- rootnode(tree_tab)$node
tree_tab[tree_tab$node == root_id,]$IGHZ <- ighz_status_root
ighz_status_mixed <- tree_tab %>% filter(IGHZ == "MIXED") %>% pull(node)
while (length(ighz_status_mixed > 0)){
  # Find and count ancestors for each node
  forebears <- lapply(ighz_status_mixed, function(n) ancestor(tree_tab, n))
  n_forebears <- sapply(forebears, nrow)
  # Pick a parent with (joint) lowest number of ancestors
  m <- sample(which(n_forebears == min(n_forebears)), 1)
  # Assign IGHZ status
  forebears_ighz <- forebears[[m]]$IGHZ
  if (any(is.na(forebears_ighz))) { ighz_status <- NA
  } else if (all(forebears_ighz == forebears_ighz[1])){ighz_status <- forebears_ighz[1]
  } else { ighz_status <- "MIXED" }
  tree_tab[tree_tab$node == ighz_status_mixed[m],]$IGHZ <- ighz_status
  # Recompute unknown node list
  ighz_status_mixed <- tree_tab %>% filter(IGHZ == "MIXED") %>% pull(node)
}

# Convert IGHZ status back to numeric
tree_tab$IGHZ <- as.numeric(tree_tab$IGHZ)

# ? Identify IGHZ-status change points and separate nodes
max_node <- max(tree_tab$node)
tree_tab <- arrange(tree_tab, node) %>% mutate(IGHZ_switch = 0)
class(tree_tab) <- c("tbl_tree", class(tree_tab))
tree_tab$IGHZ_parent <- tree_tab$IGHZ[tree_tab$parent]
ighz_status_changed <- tree_tab %>% filter(IGHZ != IGHZ_parent) %>% pull(node)
tree_tab_changed <- tree_tab %>% filter(node %in% ighz_status_changed)
tree_tab_cut <- tree_tab %>% filter(! node %in% ighz_status_changed)

# ? Insert intermediate nodes for change points
tree_tab_inter <- tree_tab_changed %>%
  mutate(node = seq(max_node+1, max_node + n()),
         IGHZ_child = IGHZ,
         IGHZ_switch = 1, IGHZ = IGHZ_parent,
         label = NA, Genus = NA, Species = NA)
tree_tab_changed$parent <- tree_tab_inter$node
tree_tab_out <- bind_rows(tree_tab_cut, tree_tab_inter, tree_tab_changed) %>%
  arrange(node)
class(tree_tab_out) <- c("tbl_tree", class(tree_tab_out))

# Convert to treedata for plotting
tree_data <- as.treedata(tree_tab_out)

#------------------------------------------------------------------------------
# Plot tree with taxon annotations
#------------------------------------------------------------------------------

tree_data@phylo$root.edge <- 0.5
tree_data@phylo$edge.length <- NA
treeplot <- suppressWarnings(ggtree(tree_data, aes(colour = IGHZ), size=1.2)) + 
  geom_rootedge(colour = colours[["CZ"]], size=1.2) +
  geom_tiplab(aes(label = paste(Genus, Species), 
                  fontface = Highlight+3), size=6.5, offset=0.3,
              family = font) +
  geom_nodepoint(aes(colour = IGHZ_child, fill = IGHZ_child, alpha = IGHZ_switch),
                 size = 5, shape = 21, stroke = 2) +
  scale_color_continuous(low = "black", high = colours[["CZ"]], lim = c(0,1),
                         labels = rev(c("Absent", "Suspected absent", "Present")),
                         breaks = sort(unique(tree_tab_out$IGHZ), decr = TRUE)) +
  scale_fill_continuous(low = "black", high = colours[["CZ"]],lim = c(0,1)) +
  scale_alpha(range = c(0, 1))
treeplot <- revts(treeplot) + suppressWarnings(xlim(NA, 6)) +
  guides(alpha = FALSE, fill = FALSE,
         colour = guide_legend(title = "IGHZ?",
                               override.aes = list(shape = 15, stroke=3))
         ) +
  theme(legend.position = "bottom",
        legend.title = element_text(family = titlefont, face = "bold",
                                    size = fontsize_base * fontscale_title),
        legend.text = element_text(family = font,
                                   size = fontsize_base  * 1.2)
  )

plot_height <- 25
plot_width <- 25

ggsave(plot = treeplot, filename = outpath, device = "svg", units = "cm",
       height=plot_height, width=plot_width)
