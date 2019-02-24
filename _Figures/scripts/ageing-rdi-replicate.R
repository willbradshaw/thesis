###############################################################################
## FIGURES & DATA                                                            ##
## Visualise cross-replicate RDI for ageing data                             ##
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
source("aux/changeo.R")
source("aux/trees.R")

# Set parameters
rdi_distance <- "euclidean" # or "cor"
rdi_iterations <- 100
rdi_constant_scale <- TRUE
rdi_transform <- TRUE
clust_method <- "average"
rootedge_length <- -0.2
treeline_width <- 1.1
tiplab_offset <- -0.1


# Configure input
# TODO: Select input settings to match other spectra in chapter
tab_path <- "../_Data/changeo/ctabs/ageing-final.tab"
seqset <- "all" # or "functional"
segments <- "VDJ" # or "VDJ"
group <- "REPLICATE" # or "INDIVIDUAL", or something else

# Output paths
filename_base <- paste0("pilot-rdi-", segments)

#------------------------------------------------------------------------------
# IMPORT FINAL TABLE (WITH COMBINED CALLS)
#------------------------------------------------------------------------------

tab <- import_tab(tab_path)

#------------------------------------------------------------------------------
# FILTER AMBIGUOUS SEGMENT CALLS
#------------------------------------------------------------------------------

has_field <- paste0("HAS_", toupper(segments))
ambig_field <- paste0(toupper(segments), "_AMBIG")

tab_filtered <- filter(tab, !!as.name(has_field), !(!!as.name(ambig_field)))

if (seqset == "functional") tab_filtered <- filter(tab_filtered, FUNCTIONAL)
if (seqset == "nonfunctional") tab_filtered <- filter(tab_filtered, !FUNCTIONAL)

#------------------------------------------------------------------------------
# COMPUTE RDI SEGMENT CALLS
#------------------------------------------------------------------------------

call_field <- paste0("BEST_", toupper(segments), "_CALL")

genes <- pull(tab, call_field)
annots <- as.character(pull(tab, group))

counts <- calcVDJcounts(genes = genes, seqAnnot = annots,
                        simplifyNames = TRUE, splitCommas = FALSE)

#------------------------------------------------------------------------------
# COMPUTE RDI DISTANCE MATRIX
#------------------------------------------------------------------------------

rdi <- calcRDI(counts, subsample = TRUE,
               distMethod = rdi_distance, nIter = rdi_iterations,
               constScale = rdi_constant_scale,
               units = ifelse(rdi_transform, "lfc", "pct"))

#------------------------------------------------------------------------------
# PERFORM CLUSTERING AND GENERATE DENDROGRAM
#------------------------------------------------------------------------------

clust <- hclust(rdi, method = clust_method)

#------------------------------------------------------------------------------
# CONVERT TO PHYLOGENETIC TREE AND VISUALISE
#------------------------------------------------------------------------------

# Convert to phylotibble and extract individual/replicate info
phylo <- as.phylo(clust)
phylo_tbl <- as_tibble(phylo) %>%
  mutate(individual = sub("(\\d-\\d\\d)(.*)", "\\1", label),
         replicate = sub("(\\d-\\d\\d)(.*)", "\\2", label),
         label = sub("(\\d-\\d\\d)(.*)", "\\1 \\2", label))
class(phylo_tbl) <- c("tbl_tree", class(phylo_tbl))

# Propagate individual info up tree
# for (i in unique(phylo_tbl %>% filter(!is.na(individual)) %>% 
#                  pull(individual))){
#   ca <- get_mrca(i, phylo)
#   of <- offspring(phylo_tbl, ca)
#   phylo_tbl[phylo_tbl$node %in% c(ca, of$node),]$individual <- i
# }
# phylo_tbl[is.na(phylo_tbl$individual),]$individual <- "GR"

# Convert back to treedata
treedata_rdi <- as.treedata(phylo_tbl)
treedata_rdi@phylo$root.edge <- rootedge_length

# Make basic plot and colour by individual
# palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
#              colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]],
#              colours[["GR"]])
theme_tre <- theme_minimal() + theme_base +
  theme(legend.position = "none",
        legend.margin = margin(b=0,t=0.5, unit="cm"),
        plot.margin = margin(t=0, unit="cm"),
        axis.title.x = element_text(margin = margin(t=0, b=0.5, unit="cm")),
        axis.title.y = element_text(margin = margin(l = 2, t = 2, unit = "cm")))


treeplot <- revts(ggtree(treedata_rdi, aes(colour = individual), 
                         size = treeline_width)) +
  geom_rootedge(colour = colours[["GR"]], size = treeline_width) +
#  scale_colour_manual(values = palette) + 
  coord_flip() + 
  scale_x_reverse(labels = function(x) -x, name = "RDI", limits = c(1.5, NA),
                  breaks = seq(0,-5)) +
  scale_y_continuous(breaks=NULL, name = "Replicate") +
  geom_tiplab(size = 3, offset = tiplab_offset, family = font,
              angle = 270) +
  theme_minimal() + theme_base + theme(legend.position = "none")

#------------------------------------------------------------------------------
# PRINCIPLE CO-ORDINATE ANALYSIS
#------------------------------------------------------------------------------

# # Perform PCoA
# pc <- pcoa(dst)
# 
# # Extract co-ordinates in first two axes into tibble for plotting
# pc_tab <- tibble(REPLICATE = rownames(as.matrix(dst)),
#                  PCO1 = pc$vectors[,1],
#                  PCO2 = pc$vectors[,2]) %>%
#   mutate(INDIVIDUAL = sub("(2-0\\d)(.*)", "\\1", REPLICATE),
#          REP = sub("(2-0\\d)(.*)", "\\2", REPLICATE))
# pc_var <- pc$values %>% pull(Broken_stick) %>% (function(x) round(x*100, 1))
# 
# g_pcoa <- ggplot(pc_tab) + geom_point(aes(x=PCO1, y=PCO2, 
#                                           colour = INDIVIDUAL), size = 3) +
#   xlab(paste0("Principle co-ordinate 1 (", pc_var[1], "%)")) +
#   ylab(paste0("Principle co-ordinate 2 (", pc_var[2], "%)")) +
#   scale_colour_manual(values = palette, name = "Individual") + 
#   theme_classic() + theme_base

#------------------------------------------------------------------------------
# COMBINE AND SAVE PLOTS
#------------------------------------------------------------------------------
# 
# plt <- plot_grid(treeplot, g_pcoa,
#                  ncol = 2, nrow = 1, labels="AUTO",
#                  label_fontfamily = titlefont, label_fontface = "plain",
#                  label_size = fontsize_base * fontscale_label)
# 
# savefig(plt, filename_base, width = 30, height = 15)
#   
