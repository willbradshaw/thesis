###############################################################################
## FIGURE                                                                    ##
## Plot diversity spectra of gut clonal repertoires                          ##
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
palette <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_group2"]],
             colours_igseq[["gut_group3"]], colours_igseq[["gut_group4"]],
             colours_igseq[["gut_group5"]])
palette_age <- c(colours_igseq[["gut_group1"]], colours_igseq[["gut_allold"]])
age_groups <- c("6", "16")
treatment_groups <- c("YI_6", "WT_16", "ABX_16", "SMT_16", "YMT_16")
qvals <- c(0,1,1.5,2,3,4)
significance_level <- 0.05

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab_grouped_group <- import_div(inpath_grouped_group)
tab_grouped_age <- import_div(inpath_grouped_age)
tab_solo_group <- import_div(inpath_solo_group)
tab_solo_age <- import_div(inpath_solo_age)

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

g_solo_age <-  ggplot(tab_solo_age) + 
  geom_line(aes(x=Q, y=D, colour = AGE_WEEKS, group = INDIVIDUAL)) + 
  geom_ribbon(aes(x=Q, ymin = D_LOWER, ymax = D_UPPER, 
                  fill = AGE_WEEKS, group = INDIVIDUAL), alpha = 0.4) +
  facet_wrap(~AGE_WEEKS, scales = "free") +
  xlab("Diversity order (q)") + 
  ylab(expression(Diversity~(""[q]*D))) +
  xlim(c(0,4)) + ylim(c(0,700)) +
  scale_colour_manual(values = palette_age, name = "Age group (weeks)") +
  scale_fill_manual(values = palette_age, name = "Age group (weeks)") +
  theme_classic() + theme_base +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm")))

g_solo_group <-  ggplot(tab_solo_group) + 
  geom_line(aes(x=Q, y=D, colour = GROUP, group = INDIVIDUAL)) + 
  geom_ribbon(aes(x=Q, ymin = D_LOWER, ymax = D_UPPER, 
                  fill = GROUP, group = INDIVIDUAL), alpha = 0.4) +
  facet_wrap(~GROUP, scales = "free") +
  xlab("Diversity order (q)") + 
  ylab(expression(Diversity~(""[q]*D))) +
  xlim(c(0,4)) + ylim(c(0,700)) +
  scale_colour_manual(values = palette, name = "Treatment group") +
  scale_fill_manual(values = palette, name = "Treatment group") +
  theme_classic() + theme_base +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm")))

g_solo <- gplot_grid(g_solo_age, g_solo_group, ncol = 1, nrow = 2)

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
                     group = !!as.name(test_by)), outlier.shape = NA) +
    geom_point(alpha = 0.4, size = 2) +
    facet_wrap(~Q, scales = "free", labeller = function(q) label_both(q, sep = " = ")) +
    scale_fill_manual(values = palette, name = x_lab) +
    ylab("Diversity") + xlab(x_lab) + ylim(c(0, NA)) +
    theme_classic() + theme_base + theme(
      legend.title = element_text(margin = margin(r = 1, unit = "cm"))
    )
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

g_solofit_both <- gplot_grid(g_solofit_age, g_solofit_group, 
                             ncol = 1, nrow = 2)
#------------------------------------------------------------------------------
# SAVE FIGURES
#------------------------------------------------------------------------------

# Alpha plots
g_alpha_both = gplot_grid(g_alpha_age, g_alpha_group)
ggsave(plot = g_solo, filename = outpath_solo, device = "svg", units = "cm",
       height=25, width = 25)
ggsave(plot = g_alpha_age, filename = outpath_alpha_age, device = "svg", units = "cm",
       height=15, width = 25)
ggsave(plot = g_alpha_group, filename = outpath_alpha_group, device = "svg", units = "cm",
       height=15, width = 25)
ggsave(plot = g_alpha_both, filename = outpath_alpha_all,  device = "svg", units = "cm",
       height = 15, width = 30)

# Boxplots
ggsave(plot = g_solofit_both, filename = outpath_solo_box,  device = "svg", units = "cm",
       height = 40, width = 40*3/4)
