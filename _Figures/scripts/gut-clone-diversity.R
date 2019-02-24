###############################################################################
## FIGURE                                                                    ##
## Plot diversity spectra of gut clonal repertoires                          ##
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

# Configure input
# TODO: Select input settings to match other spectra in chapter
seqset <- "all" # or "functional"
copy <- "NULL" # or "DUPCOUNT"
groupby <- "group" # or "age"

tab_path_grouped_group <- paste0("../_Data/changeo/spectra/gut-group",
                           "_clone-diversity-grouped_seqs-",
                           seqset, "_copy-", copy, ".tsv")
tab_path_grouped_age <- paste0("../_Data/changeo/spectra/gut-age",
                               "_clone-diversity-grouped_seqs-",
                               seqset, "_copy-", copy, ".tsv")
tab_path_solo_group <- paste0("../_Data/changeo/spectra/gut-group",
                              "_clone-diversity-solo_seqs-",
                              seqset, "_copy-", copy, ".tsv")
tab_path_solo_age <- paste0("../_Data/changeo/spectra/gut-age",
                            "_clone-diversity-solo_seqs-",
                            seqset, "_copy-", copy, ".tsv")

# Set parameters
palette <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_group2"]],
             colours_igseq[["gut_group3"]], colours_igseq[["gut_group4"]],
             colours_igseq[["gut_group5"]])
palette_age <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_allold"]])
age_groups <- c("6", "16")
treatment_groups <- c("YI_6", "WT_16", "ABX_16", "SMT_16", "YMT_16")
qvals <- c(0,1,1.5,2,3,4)
significance_level <- 0.05

# Output paths
filename_base <- "igseq-gut-clone-diversity"

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab_grouped_group <- import_div(tab_path_grouped_group)
tab_grouped_age <- import_div(tab_path_grouped_age)
tab_solo_group <- import_div(tab_path_solo_group)
tab_solo_age <- import_div(tab_path_solo_age)

# Process grouped spectra for alpha plot
tab_grouped_group <- mutate(tab_grouped_group, GROUP = sub("YI_16", "YI_6", GROUP),
                      AGE_WEEKS = factor(sub(".*_", "", GROUP), 
                                         levels = age_groups),
                      GROUP = factor(GROUP, levels = treatment_groups))
tab_grouped_age <- mutate(tab_grouped_age, 
                          AGE_WEEKS = factor(sub(".*_", "", AGE_WEEKS), 
                                          levels = age_groups))


# Process solo spectra as appropriate
tab_solo_group <- mutate(tab_solo_group, GROUP = sub("YI_16", "YI_6", GROUP),
                         AGE_WEEKS = factor(sub(".*_", "", GROUP), 
                                            levels = age_groups),
                         GROUP = factor(GROUP, levels = treatment_groups))
tab_solo_age <- mutate(tab_solo_age,
                       AGE_WEEKS = factor(sub(".*_", "", AGE_WEEKS), 
                                          levels = age_groups))

#------------------------------------------------------------------------------
# GENERATE ALPHA SPECTRA
#------------------------------------------------------------------------------

g_alpha_group <- plot_diversity_alpha(tab_grouped_group, "GROUP") +
  scale_colour_manual(values = palette, name = "Treatment group") +
  scale_fill_manual(values = palette, name = "Treatment group") +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm"))) +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
g_alpha_age <- plot_diversity_alpha(tab_grouped_age, "AGE_WEEKS") +
  scale_colour_manual(values = palette_age, name = "Age at death (weeks)") +
  scale_fill_manual(values = palette_age, name = "Age at death (weeks)") +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm"))) +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

#------------------------------------------------------------------------------
# COMPARE DIVERSITY OF DIFFERENT AGE GROUPS
#------------------------------------------------------------------------------

# Filter solo-diversity table
filter_solo_tab <- function(tab, q, test_by = "AGE_WEEKS"){
  tab %>% filter(as.character(Q) == as.character(q)) %>%
    select(!!as.name(test_by), INDIVIDUAL, D) %>%
    arrange(!!as.name(test_by), INDIVIDUAL)
}

multi_filter <- function(tab, qvals, test_by = "AGE_WEEKS"){
  bind_rows(lapply(qvals, function(q) 
    filter_solo_tab(tab, q, test_by) %>% mutate(Q = q)))
}

# Mann-Whitney U tests
mwu_diversity <- function(tab, q, test_by = "AGE_WEEKS"){
  wilcox.test(formula = as.formula(paste("D~", test_by)), 
               data = filter_solo_tab(tab, q, test_by))
}

get_mwu_grid_q_single <- function(tab, q, test_by = "AGE_WEEKS"){
  groups <- unique(levels(tab[[test_by]]))
  form <- as.formula(paste("D ~", test_by))
  L <- lapply(1:(length(groups)-1), function(a)
    lapply((a+1):length(groups), function(b)
      tab %>% filter(!!as.name(test_by) %in% c(groups[a],groups[b])) %>%
        mwu_diversity(q, test_by) %>% (function(x) x$p.value)))
  mwu_grid <- melt(L, value.name = "P") %>%
      group_by(L1) %>% 
      mutate(GROUP1 = groups[L1], GROUP2 = groups[L1 + L2]) %>%
      ungroup() %>% select(-L1, -L2) %>% mutate(Q = q)
  return(mwu_grid)
}
get_mwu_grid_q <- function(tab, qvals, test_by = "AGE_WEEKS"){
  bind_rows(lapply(qvals, function(q) 
    get_mwu_grid_q_single(tab, q, test_by)))
}

# Kruskal-Wallis non-parametric ANOVA
kruskal_diversity <- function(tab, q, test_by = "AGE_WEEKS"){
  kruskal.test(formula = as.formula(paste("D~", test_by)), 
               data = filter_solo_tab(tab, q, test_by))
}

signif_stars <- function(p){
  l <- rep("n.s.", length(p))
  l <- ifelse(p <= 0.05, "*", l)
  l <- ifelse(p <= 0.01, "**", l)
  l <- ifelse(p <= 0.001, "***", l)
  return(l)
}

# Generate plots
plot_solo_diversity <- function(tab, qvals, test_by = "AGE_WEEKS",
                                x_lab = "Age at death (weeks)", 
                                reference = age_groups, palette = palette){
  ggplot(multi_filter(tab, qvals, test_by), aes(x=!!as.name(test_by), y=D)) +
    geom_boxplot(aes(fill = factor(!!as.name(test_by), levels = reference),
                     group = !!as.name(test_by))) +
    facet_wrap(~Q, scales = "free", labeller = function(q) label_both(q, sep = " = ")) +
    scale_fill_manual(values = palette, name = x_lab) +
    ylab("Diversity") + xlab(x_lab) + ylim(c(0, NA)) +
    theme_classic() + theme_base
}

# Make boxplots and annotate with P-values
p_pos_age <- tibble(Q = qvals, XMIN = 1, XMAX = 2, XLAB = 1.5,
                Y = c(650,600,500,410,250,160), SCALE = 1)
p_age <- get_mwu_grid_q(tab_solo_age, qvals, "AGE_WEEKS") %>%
  filter(P <= significance_level) %>% mutate(LABEL = signif_stars(P)) %>%
  full_join(p_pos_age, by = "Q")
g_solofit_age <- plot_solo_diversity(tab_solo_age, qvals,
                                     palette = palette_age) +
  geom_segment(data = p_age, aes(x=XMIN, xend=XMAX, y=Y, yend = Y)) +
  geom_text(data = p_age, aes(x=XLAB, y=Y + SCALE, label = LABEL), size = 7)


#------------------------------------------------------------------------------
# COMPARE DIVERSITY OF DIFFERENT TREATMENT GROUPS
#------------------------------------------------------------------------------

p_group <- get_mwu_grid_q(tab_solo_group, qvals, "GROUP")

g_solofit_group <- plot_solo_diversity(tab_solo_group, qvals, "GROUP",
                                       "Treatment group", treatment_groups, 
                                       palette) +
  theme(axis.text.x = element_text(#size = fontsize_base * fontscale_title,
      angle = 45, hjust = 1)
  )
# No significant relationships here, so didn't add code for stars

#------------------------------------------------------------------------------
# SAVE FIGURES
#------------------------------------------------------------------------------

# Alpha plots
savefig(plot = g_alpha_age, filename = paste0(filename_base, "-alpha-age"),
        height = 15, width = 25)
savefig(plot = g_alpha_group, filename = paste0(filename_base,"-alpha-groups"),
        height = 15, width = 25)
savefig(plot = gplot_grid(g_alpha_age, g_alpha_group),#, ncol = 1, nrow = 2), 
        filename = paste0(filename_base,"-alpha"),
        height = 15, width = 30)

# Boxplots
savefig(plot = g_solofit_age, height = 20, ratio = 1.5,
        filename = paste0(filename_base, "-solo-age"))
savefig(plot = g_solofit_group, height = 20, ratio = 1.5,
        filename = paste0(filename_base, "-solo-groups"))