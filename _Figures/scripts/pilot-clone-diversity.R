###############################################################################
## FIGURE                                                                    ##
## Plot diversity spectra of pilot clonal repertoires                        ##
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
tab_path_grouped <- paste0("../_Data/changeo/spectra/",
                           "pilot_clone-diversity-grouped_seqs-",
                           seqset, "_copy-", copy, ".tsv")
tab_path_solo <- paste0("../_Data/changeo/spectra/",
                        "pilot_clone-diversity-solo_seqs-",
                        seqset, "_copy-", copy, ".tsv")
tab_path_stats_indiv <- paste0("../_Data/changeo/pilot-stats-indiv.tsv")
tab_path_stats_rep <- paste0("../_Data/changeo/pilot-stats-rep.tsv")

# Output paths
filename_base <- "pilot-clone-diversity"

# Set parameters
qvals <- c(0,1,1.5,2,3,4)
significance_level <- 0.05
plot_unit = "cm"
plot_height <- 15
plot_width <- 25

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

# Spectra
tab_grouped <- import_div(tab_path_grouped)
tab_solo <- import_div(tab_path_solo)

# Other metrics
tab_stats_indiv <- suppressMessages(read_tsv(tab_path_stats_indiv))
tab_stats_rep <- suppressMessages(read_tsv(tab_path_stats_rep))

#------------------------------------------------------------------------------
# GENERATE ALPHA AND BETA SPECTRA
#------------------------------------------------------------------------------

palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
             colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]])

g_alpha <- plot_diversity_alpha(tab_grouped, "INDIVIDUAL") +
  scale_colour_manual(values = palette, name = "Individual") +
  scale_fill_manual(values = palette, name = "Individual")
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
  xlab("Diversity order (q)") + 
  ylab(expression(Diversity~(""[q]*D))) +
  theme_classic() + theme_base


#------------------------------------------------------------------------------
# COMPARE DIVERSITY OF DIFFERENT AGE GROUPS
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
           get_mwu_grid(test_by, test_on) %>%
           process_mwu_grid(ymax, reference, signif_level, scale) %>%
           mutate(Q = q)) %>% bind_rows
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
                     group = !!as.name(test_by))) +
    facet_wrap(~Q, scales = "free", labeller = function(q) label_both(q, sep = " = ")) +
    scale_fill_manual(values = palette, name = x_lab) +
    ylab("Diversity") + xlab(x_lab) + ylim(c(0, NA)) +
    theme_classic() + theme_base
}

# Make boxplots and annotate with P-values
kruskal_yvals <- multi_kruskal(tab_solo, qvals, "INDIVIDUAL", TRUE) %>%
  mutate(Y = c(100, 80, 50, 20, 10, 5), 
         LABEL = paste0("P(KWT) = ", signif(P,3)))
mwu <- multi_mwu_grid(tab_solo, qvals, "INDIVIDUAL", "D", 2000,
                      tab_solo %>% pull(INDIVIDUAL) %>% unique,
                      significance_level, 0.2)
g_indiv_alpha <- plot_solo_diversity(tab_solo, qvals, "INDIVIDUAL", "Individual",
                                     tab_solo %>% pull(INDIVIDUAL) %>% unique,
                                     palette) +
  geom_text(aes(y=Y, label=LABEL), data = kruskal_yvals, x=0.8, size = 3.8, 
            family = font, hjust = 0)
# Not adding MWU data since nothing is significant

#------------------------------------------------------------------------------
# RELATE
#------------------------------------------------------------------------------

# Compute correlations between D and other metrics at each Q-value
rtab_indiv <- tab_grouped %>% 
  full_join(tab_stats_indiv, by = "INDIVIDUAL") %>% group_by(Q) %>% 
  summarise(R_P20 = cor(D, P20_Observed), R_S_Filtered = cor(D, S_Filtered), 
            R_S_All = cor(D, S_All))
rtab_rep <- tab_solo %>% 
  full_join(tab_stats_rep, by = "REPLICATE") %>% group_by(Q) %>% 
  summarise(R_P20 = cor(D, P20_Observed), R_S_Filtered = cor(D, S_Filtered), 
            R_S_All = cor(D, S_All))

# Melt for plotting
rtab_indiv_melt <- rtab_indiv %>% melt(id.vars = "Q", value.name = "r") %>%
  mutate(variable = sub("^R_", "", variable))
rtab_rep_melt <- rtab_rep %>% melt(id.vars = "Q", value.name = "r") %>%
  mutate(variable = sub("^R_", "", variable))

# Make plots
rplot_indiv <- ggplot(rtab_indiv_melt) + 
  geom_line(aes(x=Q, y=r, colour = variable)) +
  xlab("Diversity order (q)") + ylim(c(-1,0)) +
  ylab(expression("Pearson correlation with"~scriptstyle(""[q]*D^alpha))) +
  scale_colour_discrete(name = "Metric",
                        labels = c("P20", "Zipf exponent (all)", 
                                   "Zipf exponent (filtered)")) +
  theme_classic() + theme_base
rplot_rep <- ggplot(rtab_rep_melt) + 
  geom_line(aes(x=Q, y=r, colour = variable)) +
  xlab("Diversity order (q)") + ylim(c(-1,0)) + 
  ylab(expression("Pearson correlation with"~""[q]*D)) +
  scale_colour_discrete(name = "Metric",
                        labels = c("P20", "Zipf exponent (all)", 
                                   "Zipf exponent (filtered)")) +
  theme_classic() + theme_base


#------------------------------------------------------------------------------
# SAVE FIGURES
#------------------------------------------------------------------------------

savefig(spectra, filename_base, height = plot_height, width = plot_width)
savefig(g_solo, paste0(filename_base, "-solo-spectra"),
        height = 15, ratio = 1.3)
savefig(g_indiv_alpha, paste0(filename_base, "-solo-box"),
        height = 20, ratio = 1.5)
savefig(rplot_rep, paste0(filename_base, "-metrics-cor"),
        height = 15, ratio = 1.3)