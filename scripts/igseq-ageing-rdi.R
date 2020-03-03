###############################################################################
## FIGURES & DATA                                                            ##
## Visualise cross-individual RDI for ageing data                            ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Set parameters
rdi_distance <- "euclidean" # or "cor"
rdi_iterations <- 100
rdi_constant_scale <- TRUE
rdi_transform <- TRUE
clust_method <- "average"
rootedge_length <- -0.2
treeline_width <- 1.1
tiplab_offset <- -0.1
palette <- c(colours_igseq[["ageing_group1"]], colours_igseq[["ageing_group2"]],
             colours_igseq[["ageing_group3"]], colours_igseq[["ageing_group4"]])
age_groups <- c("39", "56", "73", "128")


# Configure input
seqset <- "all" # or "functional"
segments <- "VJ" # or "VDJ"
group <- "INDIVIDUAL"

#------------------------------------------------------------------------------
# IMPORT FINAL TABLE (WITH COMBINED CALLS)
#------------------------------------------------------------------------------

tab <- import_tab(inpath)

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

genes <- pull(tab_filtered, call_field)
annots <- as.character(pull(tab_filtered, group))

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
# PRINCIPLE CO-ORDINATE ANALYSIS
#------------------------------------------------------------------------------

indivs <- tab %>% group_by(AGE_DAYS, INDIVIDUAL) %>% summarise %>% ungroup

# Perform PCoA
pc <- pcoa(rdi)

# Extract co-ordinates in first two axes into tibble for plotting
pc_tab <- tibble(REPLICATE = rownames(as.matrix(rdi)),
                 PCO1 = pc$vectors[,1],
                 PCO2 = pc$vectors[,2]) %>%
  mutate(INDIVIDUAL = sub("(2-0\\d)(.*)", "\\1", REPLICATE),
         REP = sub("(2-0\\d)(.*)", "\\2", REPLICATE)) %>%
  full_join(indivs, by = "INDIVIDUAL")

pc_var <- pc$values %>% pull(Broken_stick) %>% (function(x) round(x*100, 1))

g_pcoa_all <- ggplot(pc_tab) + 
  geom_point(aes(x=PCO1, y=PCO2, colour = factor(AGE_DAYS, levels = age_groups))) +
  xlab(paste0("Principle co-ordinate 1 (", pc_var[1], "%)")) +
  ylab(paste0("Principle co-ordinate 2 (", pc_var[2], "%)")) +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_x_continuous(limits=c(-2.6,2.8)) + scale_y_continuous(limits=c(-2,2.6)) +
  theme_classic() + theme_base +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm")))
g_pcoa_facet <- g_pcoa_all +
  facet_wrap(~factor(AGE_DAYS, levels = age_groups), scales = "free")


#------------------------------------------------------------------------------
# VISUALISE INTRA-GROUP DISTANCES
#------------------------------------------------------------------------------

indivs <- tab %>% group_by(AGE_DAYS, INDIVIDUAL) %>% summarise %>%
  ungroup() %>% arrange(INDIVIDUAL) %>% mutate(N = row_number())
indiv_index <- setNames(indivs$N, indivs$INDIVIDUAL)
indiv_ages  <-  setNames(indivs$AGE_DAYS, indivs$INDIVIDUAL)

dist_tab <- rdi %>% as.matrix %>% 
  melt(varnames = c("ID1", "ID2"), value.name = "RDI") %>% as.tibble %>%
  mutate(ROW = indiv_index[ID1], COL = indiv_index[ID2],
         AGE1 = indiv_ages[ID1], AGE2 = indiv_ages[ID2])

# Remove redundant distances and restrict to matching age groups
dist_tab_filtered <- filter(dist_tab, ROW > COL, AGE1 == AGE2) %>%
  select(AGE_DAYS = AGE1, ID1, ID2, RDI)

# Get nearest-neighbour distances for each individual
dist_tab_nn <- dist_tab %>% filter(ID1 != ID2, AGE1 == AGE2) %>%
  group_by(ID1) %>% filter(RDI == min(RDI)) %>% 
  select(AGE_DAYS = AGE1, INDIVIDUAL = ID1, NEIGHBOUR = ID2, RDI)

# Make plots
g_intra_all <- ggplot(dist_tab_filtered, aes(x=as.numeric(AGE_DAYS), y=RDI)) +
  geom_boxplot(aes(fill = factor(AGE_DAYS, levels = age_groups)), 
               outlier.shape = NA) +
  geom_point(size = 2, alpha = 0.3) +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  xlab("Age at death (days)") + 
  theme_classic() + theme_base + 
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm")))
g_intra_nn <- ggplot(dist_tab_nn, aes(x=as.numeric(AGE_DAYS), y=RDI)) +
  geom_boxplot(aes(fill = factor(AGE_DAYS, levels = age_groups)), 
               outlier.shape = NA) +
  geom_point(size = 2, alpha = 0.4) +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  xlab("Age at death (days)") + ylab("Nearest-neighbour RDI") +
  theme_classic() + theme_base +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm")))

#------------------------------------------------------------------------------
# TEST FOR AGE EFFECT AND ANNOTATE PLOTS
#------------------------------------------------------------------------------

# Kruskal-Wallis tests
kruskal_test_all <- kruskal.test(RDI ~ as.numeric(AGE_DAYS), dist_tab_filtered)
kruskal_test_nn <- kruskal.test(RDI ~ as.numeric(AGE_DAYS), dist_tab_nn)

# Mann-Whitney-U tests
dist_max_all <- dist_tab_filtered %>% pull(RDI) %>% max
dist_max_nn <- dist_tab_nn %>% pull(RDI) %>% max

signif_stars <- function(p){
  l <- rep("n.s.", length(p))
  l <- ifelse(p <= 0.05, "*", l)
  l <- ifelse(p <= 0.01, "**", l)
  l <- ifelse(p <= 0.001, "***", l)
  return(l)
}
pull_rdis <- function(dtab, age_group) dtab %>% 
  filter(AGE_DAYS == age_group) %>% pull(RDI)
get_mwu_grid <- function(dtab){
  age_groups <- unique(dtab[["AGE_DAYS"]])
  L <- lapply(1:(length(age_groups)-1), function(a) 
    lapply((a+1):length(age_groups), function(b) wilcox.test(
      pull_rdis(dtab, age_groups[a]), pull_rdis(dtab, age_groups[b]))$p.value))
  mwu_grid <- melt(L, value.name = "P") %>%
    group_by(L1) %>% 
    mutate(AGE1 = age_groups[L1], AGE2 = age_groups[L1 + L2]) %>%
    ungroup() %>% select(-L1, -L2)
  return(mwu_grid)
}

process_mwu_grid <- function(mwu_grid, dist_max, signif_level = 0.05, 
                                 scale = 0.2){
  mwu_grid %>% filter(P <= signif_level) %>%
    mutate(AGE1 = as.numeric(AGE1), AGE2 = as.numeric(AGE2),
           AGE_AVG = (AGE1 + AGE2)/2, AGE_MIN = pmin(AGE1,AGE2),
           AGE_DIFF = pmax(AGE1, AGE2) - pmin(AGE1,AGE2)) %>%
    arrange(AGE_MIN, AGE_DIFF) %>%
    mutate(Y_BAR = dist_max + scale * (row_number()-0.5),
           Y_LAB = Y_BAR + scale/4,
           LABEL = signif_stars(P))
}

# Prepare tibbles of significance annotations
mwu_grid_all <- dist_tab_filtered %>% get_mwu_grid %>% 
  process_mwu_grid(dist_max_all)
mwu_grid_nn <- dist_tab_nn %>% get_mwu_grid %>% 
  process_mwu_grid(dist_max_nn, scale = 0.12)
kwt_grid_all <- tibble(X = 120, Y = 3.6, 
                       P = signif(kruskal_test_all$p.value, 3),
                       LAB = paste0("P(KWT) = ", P))
kwt_grid_nn <- tibble(X = 120, Y = 3.6, 
                       P = signif(kruskal_test_nn$p.value, 3),
                       LAB = paste0("P(KWT) = ", format(P, scientific=TRUE)))

# Annotate plots with test results
g_intra_all_annot <- g_intra_all +
  geom_segment(aes(x=AGE1,xend=AGE2,y=Y_BAR,yend=Y_BAR), data = mwu_grid_all) + 
  geom_text(aes(x=AGE_AVG, y=Y_LAB, label=LABEL), data = mwu_grid_all, size=6) +
  geom_text(aes(x=X,y=Y,label=LAB), family=font, size=3.5, data = kwt_grid_all)
g_intra_nn_annot <- g_intra_nn +
  geom_segment(aes(x=AGE1,xend=AGE2,y=Y_BAR,yend=Y_BAR), data = mwu_grid_nn) + 
  geom_text(aes(x=AGE_AVG, y=Y_LAB, label=LABEL), data = mwu_grid_nn, size=6) +
  geom_text(aes(x=X,y=Y,label=LAB), family=font, size=3.5, data = kwt_grid_nn)

#------------------------------------------------------------------------------
# COMBINE FIGURES WITH SINGLE LEGEND
#------------------------------------------------------------------------------

# Extract legend
g <- ggplotGrob(g_pcoa_all + 
                  guides(colour = guide_legend(override.aes = list(size=4))))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)

# Combine plots without legend
plt_intra <- plot_grid(g_intra_all_annot + theme(legend.position = "none"),
                 g_intra_nn_annot + theme(legend.position = "none"),
                 g_pcoa_all + theme(legend.position = "none"),
                 g_pcoa_facet + theme(legend.position = "none"),
                 ncol = 2, nrow = 2, labels="AUTO",
                 label_fontfamily = titlefont, label_fontface = "plain",
                 label_size = fontsize_base * fontscale_label)
combined <- arrangeGrob(plt_intra,
                        legend,
                        ncol = 1, nrow = 2,
                        heights = unit.c(unit(1, "npc") - lheight, lheight))

# Visualise plot
plot_unit = "cm"
plot_height <- 27
plot_width <- 25
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
grid.draw(combined)
popViewport(1)

plt_out <- grid.grab()

#------------------------------------------------------------------------------
# COMBINE AND SAVE PLOTS
#------------------------------------------------------------------------------

ggsave(plot = plt_out, filename = outpath, device = "svg", units = "cm",
       height=plot_height, width = plot_width)
