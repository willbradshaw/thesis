###############################################################################
## FIGURE                                                                    ##
## Plot diversity spectra of pilot VJ repertoires                            ##
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
palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
             colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]])
palette_ext <- c(palette, colours[["X1"]], colours[["X2"]], colours[["LR"]])
individuals <- paste0("2-0", seq(3,6))
qvals <- c(0,1,1.5,2,3,4)
significance_level <- 0.05
plot_unit = "cm"
plot_height <- 15
plot_width <- 25

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

# Spectra
tab_grouped <- import_div(inpath_grouped)
tab_solo <- import_div(inpath_solo)

#------------------------------------------------------------------------------
# GENERATE ALPHA AND BETA SPECTRA
#------------------------------------------------------------------------------

g_alpha <- plot_diversity_alpha(tab_grouped, "INDIVIDUAL") +
  scale_colour_manual(values = palette, name = "Individual") +
  scale_fill_manual(values = palette, name = "Individual") +
  theme(legend.title = element_text(margin=margin(r=0.5, unit="cm")))
g_beta <- plot_diversity_beta_scaled(tab_grouped, "INDIVIDUAL") +
  scale_colour_manual(values = palette, name = "Individual") +
  scale_fill_manual(values = palette, name = "Individual")

spectra <- gplot_grid_onelegend(g_alpha, g_beta, plot_height = plot_height,
                                plot_width = plot_width, plot_unit = plot_unit,
                                ncol = 2, nrow = 1)

g_solo <-  ggplot(tab_solo) + 
  geom_line(aes(x=Q, y=D, colour = INDIVIDUAL, group = REPLICATE)) + 
  geom_ribbon(aes(x=Q, ymin = D_LOWER, ymax = D_UPPER, 
                  fill = INDIVIDUAL, group = REPLICATE), alpha = 0.4) +
  facet_wrap(~INDIVIDUAL, scales = "free") +
  xlab("Diversity order (q)") + 
  ylab(expression(Diversity~(""[q]*D))) +
  xlim(c(0,4)) + ylim(c(0,150)) +
  scale_colour_manual(values = palette, name = "Individual") +
  scale_fill_manual(values = palette, name = "Individual") +
  theme_classic() + theme_base

#------------------------------------------------------------------------------
# COMPARE DIVERSITY OF DIFFERENT INDIVIDUALS
#------------------------------------------------------------------------------

# Filter solo-diversity table
filter_solo_tab <- function(tab, q, test_by = "INDIVIDUAL"){
  tab %>% filter(as.character(Q) == as.character(q)) %>%
    select(!!as.name(test_by), INDIVIDUAL, D) %>%
    arrange(!!as.name(test_by), INDIVIDUAL)
}

multi_filter <- function(tab, qvals, test_by = "INDIVIDUAL"){
  bind_rows(lapply(qvals, function(q) 
    filter_solo_tab(tab, q, test_by) %>% mutate(Q = q)))
}

# Mann-Whitney-U tests
signif_stars <- function(p){
  l <- rep("n.s.", length(p))
  l <- ifelse(p <= 0.05, "*", l)
  l <- ifelse(p <= 0.01, "**", l)
  l <- ifelse(p <= 0.001, "***", l)
  return(l)
}

get_mwu_grid <- function(dtab, test_by = "INDIVIDUAL", test_on = "D"){
  groups <- unique(dtab[[test_by]])
  L <- lapply(1:(length(groups)-1), function(a) 
    lapply((a+1):length(groups), function(b) wilcox.test(
      dtab %>% filter(!!as.name(test_by) == groups[a]) %>% pull(test_on),
      dtab %>% filter(!!as.name(test_by) == groups[b]) %>% pull(test_on))$p.value))
  mwu_grid <- melt(L, value.name = "P") %>%
    group_by(L1) %>% 
    mutate(GRP1 = groups[L1], GRP2 = groups[L1 + L2]) %>%
    ungroup() %>% select(-L1, -L2)
  return(mwu_grid)
}

signif_stars <- function(p){
  l <- rep("n.s.", length(p))
  l <- ifelse(p <= 0.05, "*", l)
  l <- ifelse(p <= 0.01, "**", l)
  l <- ifelse(p <= 0.001, "***", l)
  return(l)
}

process_mwu_grid <- function(mwu_grid, ymax, reference,
                             signif_level = 0.05, scale = 0.2){
  mwu_grid %>% filter(P <= signif_level) %>%
    mutate(POS1 = match(GRP1, reference), POS2 = match(GRP2, reference),
           POS_AVG = (POS1 + POS2)/2, POS_MIN = pmin(POS1, POS2),
           POS_DIFF = pmax(POS1,POS2) - pmin(POS1,POS2)) %>%
    arrange(POS_MIN, POS_DIFF) %>% # group by something here?
    mutate(Y_BAR = ymax + scale * (row_number()-0.5),
           Y_LAB = Y_BAR + scale/4,
           LABEL = signif_stars(P))
}

multi_mwu_grid <- function(tab, qvals, test_by = "INDIVIDUAL", test_on = "D",
                           ymax, reference, signif_level = 0.05, 
                           scale = 0.2){
  lapply(qvals, function(q) tab %>% filter_solo_tab(q, test_by) %>%
           get_mwu_grid(test_by, test_on) %>% mutate(Q = q)) %>% bind_rows %>%
    process_mwu_grid(ymax, reference, signif_level, scale)
} # TODO: Allow ymax to be set by Q-value

# Kruskal-Wallis non-parametric ANOVA
kruskal_diversity <- function(tab, q, test_by = "INDIVIDUAL",
                              as_factor = FALSE){
  ifelse(as_factor, paste0("as.factor(", test_by, ")"), test_by) %>% 
    paste0("D~", .) %>% as.formula(.) %>%
    kruskal.test(formula = ., data = filter_solo_tab(tab, q, test_by))
}
multi_kruskal <- function(tab, qvals, test_by = "INDIVIDUAL",
                          as_factor = FALSE){
  tibble(Q = qvals, P = sapply(qvals, function(q) 
    kruskal_diversity(tab, q, test_by, as_factor)$p.value))
}

# Generate plots
plot_solo_diversity <- function(tab, qvals, test_by = "INDIVIDUAL",
                                x_lab = "INDIVIDUAL", 
                                reference = age_groups, palette = palette){
  ggplot(multi_filter(tab, qvals, test_by), aes(x=!!as.name(test_by), y=D)) +
    geom_boxplot(aes(fill = factor(!!as.name(test_by), levels = reference),
                     group = !!as.name(test_by)), outlier.shape = NA) +
    geom_point(alpha = 0.4, size = 2) +
    facet_wrap(~Q, scales = "free", labeller = function(q) label_both(q, sep = " = ")) +
    scale_fill_manual(values = palette, name = x_lab) +
    ylab("Diversity") + xlab(x_lab) + ylim(c(0, NA)) +
    theme_classic() + theme_base
}

# Make boxplots and annotate with P-values
kruskal_yvals <- multi_kruskal(tab_solo, qvals, "INDIVIDUAL", TRUE) %>%
  mutate(Y = c(15, 7, 5, 5, 4, 3.5), 
         LABEL = paste0("P(KWT) = ", signif(P,3)))
mwu <- multi_mwu_grid(tab_solo, qvals, "INDIVIDUAL", "D", 2000, individuals,
                      significance_level, 0.2)
g_indiv_alpha <- plot_solo_diversity(tab_solo, qvals, "INDIVIDUAL", "Individual",
                                     tab_solo %>% pull(INDIVIDUAL) %>% unique,
                                     palette) +
  geom_text(aes(y=Y, label=LABEL), data = kruskal_yvals, x=0.8, size = 3.8, 
            family = font, hjust = 0)

#------------------------------------------------------------------------------
# SAVE FIGURES
#------------------------------------------------------------------------------

ggsave(plot = spectra, filename = outpath_spectra_grouped, device = "svg", units = "cm",
       height=plot_height, width = plot_width)
ggsave(plot = g_solo, filename = outpath_spectra_solo, device = "svg", units = "cm",
       height=15, width = 20)
ggsave(plot = g_indiv_alpha, filename = outpath_box_solo, device = "svg", units = "cm",
       height=20, width = 30)
