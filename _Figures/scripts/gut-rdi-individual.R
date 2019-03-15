###############################################################################
## FIGURES & DATA                                                            ##
## Visualise cross-individual RDI for ageing data                            ##
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
palette <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_group2"]],
             colours_igseq[["gut_group3"]], colours_igseq[["gut_group4"]],
             colours_igseq[["gut_group5"]])
palette_age <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_allold"]])
age_groups <- c("6", "16")
treatment_groups <- c("YI_6", "WT_16", "ABX_16", "SMT_16", "YMT_16")
individuals_excluded <- c("1274", "1309")
plot_height <- 15
plot_width <- 25

# Configure input
# TODO: Select input settings to match other spectra in chapter
tab_path <- "../_Data/changeo/ctabs/gut-final.tab"
seqset <- "all" # or "functional"
segments <- "VJ" # or "VDJ"
group <- "INDIVIDUAL"

# Output paths
filename_base <- paste0("igseq-gut-rdi-", segments, "-", tolower(group))

#------------------------------------------------------------------------------
# IMPORT FINAL TABLE (WITH COMBINED CALLS)
#------------------------------------------------------------------------------

tab <- import_tab(tab_path) %>% 
  mutate(GROUP = sub("YI_16", "YI_6", GROUP))  %>%
  filter(!INDIVIDUAL %in% individuals_excluded)


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

indivs <- tab %>% group_by(AGE_WEEKS, GROUP, INDIVIDUAL) %>% 
  summarise %>% ungroup

# Perform PCoA
pc <- pcoa(rdi)

# Extract co-ordinates in first two axes into tibble for plotting
pc_tab <- tibble(INDIVIDUAL = rownames(as.matrix(rdi)),
                 PCO1 = pc$vectors[,1],
                 PCO2 = pc$vectors[,2]) %>%
  mutate(INDIVIDUAL = sub("(2-0\\d)(.*)", "\\1", INDIVIDUAL),
         REP = sub("(2-0\\d)(.*)", "\\2", INDIVIDUAL)) %>%
  full_join(indivs, by = "INDIVIDUAL")

pc_var <- pc$values %>% pull(Broken_stick) %>% (function(x) round(x*100, 1))
xrange <- c(floor(min(pc_tab$PCO1) * 5), ceiling(max(pc_tab$PCO1) * 5)) / 5
yrange <- c(floor(min(pc_tab$PCO2) * 5), ceiling(max(pc_tab$PCO2) * 5)) / 5

# Plot on single set of axes
g_pcoa_base <- ggplot(pc_tab, aes(x=PCO1, y=PCO2)) +
  xlab(paste0("Principle co-ordinate 1 (", pc_var[1], "%)")) +
  ylab(paste0("Principle co-ordinate 2 (", pc_var[2], "%)")) +
  scale_x_continuous(limits=xrange) + scale_y_continuous(limits=yrange) +
  theme_classic() + theme_base +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm")))
g_pcoa_age_all <- g_pcoa_base +
  geom_point(aes(colour = factor(AGE_WEEKS, levels = age_groups)), size = 3) +
  scale_colour_manual(values = palette_age, name = "Age group (weeks)")
g_pcoa_group_all <- g_pcoa_base +
  geom_point(aes(colour = factor(GROUP, levels = treatment_groups)), size = 3) +
  scale_colour_manual(values = palette, name = "Treatment group")
  
# Facet PCOA plots
g_pcoa_age_facet <- g_pcoa_age_all +
  facet_wrap(~factor(AGE_WEEKS, levels = age_groups), scales = "free")
g_pcoa_group_facet <- g_pcoa_group_all +
  facet_wrap(~factor(GROUP, levels = treatment_groups), scales = "free")

#------------------------------------------------------------------------------
# VISUALISE INTRA-GROUP DISTANCES
#------------------------------------------------------------------------------

indivs <- tab %>% group_by(GROUP, AGE_WEEKS, INDIVIDUAL) %>% summarise %>%
  ungroup() %>% arrange(INDIVIDUAL) %>% mutate(N = row_number())
indiv_index <- setNames(indivs$N, indivs$INDIVIDUAL)
indiv_ages  <-  setNames(indivs$AGE_WEEKS, indivs$INDIVIDUAL)
indiv_groups <- setNames(indivs$GROUP, indivs$INDIVIDUAL)

dist_tab <- rdi %>% as.matrix %>% 
  melt(varnames = c("ID1", "ID2"), value.name = "RDI") %>% as.tibble %>%
  mutate(ROW = indiv_index[ID1], COL = indiv_index[ID2],
         AGE1 = indiv_ages[ID1], AGE2 = indiv_ages[ID2],
         GRP1 = indiv_groups[ID1], GRP2 = indiv_groups[ID2])

# Remove redundant distances and restrict to matching age or treatment groups
dist_tab_filtered_age <- filter(dist_tab, ROW > COL, AGE1 == AGE2) %>%
  select(AGE_WEEKS = AGE1, ID1, ID2, RDI)
dist_tab_filtered_group <- filter(dist_tab, ROW > COL, GRP1 == GRP2) %>%
  select(GROUP = GRP1, ID1, ID2, RDI)

# Get nearest-neighbour distances for each individual
dist_tab_nn_age <- dist_tab %>% filter(ID1 != ID2, AGE1 == AGE2) %>%
  group_by(ID1) %>% filter(RDI == min(RDI)) %>% 
  select(AGE_WEEKS = AGE1, INDIVIDUAL = ID1, NEIGHBOUR = ID2, RDI)
dist_tab_nn_group <- dist_tab %>% filter(ID1 != ID2, GRP1 == GRP2) %>%
  group_by(ID1) %>% filter(RDI == min(RDI)) %>% 
  select(GROUP = GRP1, INDIVIDUAL = ID1, NEIGHBOUR = ID2, RDI)


# Make plots
plot_rdi_intra <- function(dist_tab, test_by = "AGE_WEEKS",
                           reference = age_groups, palette = palette_age,
                           group_name = "Age group (weeks)",
                           point_size = 3, point_alpha = 0.4){
  ggplot(dist_tab, aes(x=factor(!!as.name(test_by), levels = reference), y=RDI)) +
    geom_boxplot(aes(fill = factor(!!as.name(test_by), levels = reference)), 
                 outlier.shape = NA) +
    geom_point(size = point_size, alpha = point_alpha) +
    xlab(group_name) +
    scale_colour_manual(values = palette, name = group_name) +
    scale_fill_manual(values = palette, name = group_name) +
    theme_classic() + theme_base + theme(legend.position = "none")
}

g_intra_age_all <- plot_rdi_intra(dist_tab_filtered_age, point_size = 2,
                                  point_alpha = 0.3)
g_intra_age_nn <- plot_rdi_intra(dist_tab_nn_age, point_size = 2, 
                                 point_alpha = 0.3) +
  ylab("Nearest-neighbour RDI")
g_intra_group_all <- plot_rdi_intra(dist_tab_filtered_group, "GROUP",
                                    treatment_groups, palette,
                                    "Treatment group")
g_intra_group_nn <- plot_rdi_intra(dist_tab_nn_group, "GROUP",
                                    treatment_groups, palette, 
                                   "Treatment group") +
  ylab("Nearest-neighbour RDI")

#------------------------------------------------------------------------------
# TEST FOR AGE/GROUP EFFECT AND ANNOTATE PLOTS
#------------------------------------------------------------------------------

# Mann-Whitney-U tests
dist_max_all_age <- dist_tab_filtered_age %>% pull(RDI) %>% max
dist_max_nn_age <- dist_tab_nn_age %>% pull(RDI) %>% max
dist_max_all_group <- dist_tab_filtered_group %>% pull(RDI) %>% max
dist_max_nn_group <- dist_tab_nn_group %>% pull(RDI) %>% max

signif_stars <- function(p){
  l <- rep("n.s.", length(p))
  l <- ifelse(p <= 0.05, "*", l)
  l <- ifelse(p <= 0.01, "**", l)
  l <- ifelse(p <= 0.001, "***", l)
  return(l)
}
pull_rdis <- function(dtab, group, test_by) dtab %>% 
  filter(!!as.name(test_by) == group) %>% pull(RDI)

get_mwu_grid <- function(dtab, test_by){
  groups <- unique(dtab[[test_by]])
  L <- lapply(1:(length(groups)-1), function(a) 
    lapply((a+1):length(groups), function(b) wilcox.test(
      pull_rdis(dtab, groups[a], test_by), 
      pull_rdis(dtab, groups[b], test_by))$p.value))
  mwu_grid <- melt(L, value.name = "P") %>%
    group_by(L1) %>% 
    mutate(GRP1 = groups[L1], GRP2 = groups[L1 + L2]) %>%
    ungroup() %>% select(-L1, -L2)
  return(mwu_grid)
}

process_mwu_grid <- function(mwu_grid, dist_max, reference = age_groups,
                             signif_level = 0.05, scale = 0.2){
  mwu_grid %>% filter(P <= signif_level) %>%
    mutate(POS1 = match(GRP1, reference), POS2 = match(GRP2, reference),
           POS_AVG = (POS1 + POS2)/2, POS_MIN = pmin(POS1, POS2),
           POS_DIFF = pmax(POS1,POS2) - pmin(POS1,POS2)) %>%
    arrange(POS_MIN, POS_DIFF) %>% # group by something here?
    mutate(Y_BAR = dist_max + scale * (row_number()-0.5),
           Y_LAB = Y_BAR + scale/4,
           LABEL = signif_stars(P))
}

# Prepare tibbles of significance annotations
mwu_grid_all_age <- dist_tab_filtered_age %>% get_mwu_grid("AGE_WEEKS") %>% 
  process_mwu_grid(dist_max_all_age, age_groups, scale = 0.5)
mwu_grid_all_group <- dist_tab_filtered_group %>% get_mwu_grid("GROUP") %>% 
  process_mwu_grid(dist_max_all_group, treatment_groups, scale = 0.5)
mwu_grid_nn_age <- dist_tab_nn_age %>% get_mwu_grid("AGE_WEEKS") %>% 
  process_mwu_grid(dist_max_nn_age, age_groups, scale = 0.4)
mwu_grid_nn_group <- dist_tab_nn_group %>% get_mwu_grid("GROUP") %>% 
  process_mwu_grid(dist_max_nn_group, treatment_groups, scale = 0.4)

# Annotate plots with test results
g_intra_age_all_annot <- g_intra_age_all +
  geom_segment(aes(x=POS1,xend=POS2,y=Y_BAR,yend=Y_BAR), data = mwu_grid_all_age) + 
  geom_text(aes(x=POS_AVG, y=Y_LAB, label=LABEL), data = mwu_grid_all_age, size=8)
g_intra_group_all_annot <- g_intra_group_all +
  geom_segment(aes(x=POS1,xend=POS2,y=Y_BAR,yend=Y_BAR), data = mwu_grid_all_group) + 
  geom_text(aes(x=POS_AVG, y=Y_LAB, label=LABEL), data = mwu_grid_all_group, size=8)
g_intra_age_nn_annot <- g_intra_age_nn +
  geom_segment(aes(x=POS1,xend=POS2,y=Y_BAR,yend=Y_BAR), data = mwu_grid_nn_age) + 
  geom_text(aes(x=POS_AVG, y=Y_LAB, label=LABEL), data = mwu_grid_nn_age, size=8)
g_intra_group_nn_annot <- g_intra_group_nn +
  geom_segment(aes(x=POS1,xend=POS2,y=Y_BAR,yend=Y_BAR), data = mwu_grid_nn_group) + 
  geom_text(aes(x=POS_AVG, y=Y_LAB, label=LABEL), data = mwu_grid_nn_group, size=8)

#------------------------------------------------------------------------------
# COMBINE FIGURE GROUPS WITH SINGLE LEGEND AND SAVE
#------------------------------------------------------------------------------

g_age_all <- gplot_grid_onelegend(g_pcoa_age_all, g_pcoa_age_facet, 
                                  g_intra_age_all_annot, g_intra_age_nn_annot, 
                                  plot_height = plot_height * 2, 
                                  plot_width = plot_width, nrow = 2)
g_group_all <- gplot_grid_onelegend(g_pcoa_group_all, g_pcoa_group_facet, 
                                  g_intra_group_all_annot, g_intra_group_nn_annot, 
                                  plot_height = plot_height * 2, 
                                  plot_width = plot_width, nrow = 2)

savefig(g_age_all, paste0(filename_base, "-age"),
        height = plot_height * 2, width = plot_width)
savefig(g_group_all, paste0(filename_base, "-group"),
        height = plot_height * 2, width = plot_width)