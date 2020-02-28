###############################################################################
## FIGURE                                                                    ##
## Plot diversity spectra of pilot clonal repertoires                        ##
###############################################################################
aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

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
tab_grouped <- import_div(inpath_grouped)
tab_solo <- import_div(inpath_solo)

# Other metrics
tab_stats_indiv <- suppressMessages(read_tsv(inpath_stats_indiv))
tab_stats_rep <- suppressMessages(read_tsv(inpath_stats_rep))

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
  facet_wrap(~INDIVIDUAL, scales = "free") +
  xlab("Diversity order (q)") + 
  ylab(expression(Diversity~(""[q]*D))) +
  xlim(c(0,4)) + ylim(c(0,2000)) +
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
    geom_point(alpha = 0.4, size = 2) +
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
# RELATE DIVERSITY TO EXPANSION METRICS
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

# Compute individual correlations from average over replicates
rtab_rep_avg <- tab_solo %>% full_join(tab_stats_rep, by = "REPLICATE") %>%
  group_by(Q, INDIVIDUAL) %>%
  summarise(D = ifelse(first(Q) != 1,
                       mean(D^(1-first(Q)))^(1/(1-first(Q))),
                       exp(mean(log(D)))),
            S_Filtered = mean(S_Filtered), S_All = mean(S_All), 
            P20_Observed = mean(P20_Observed)) %>%
  summarise(R_P20 = cor(D, P20_Observed), R_S_Filtered = cor(D, S_Filtered), 
            R_S_All = cor(D, S_All))

# Melt for plotting
melt_rtab <- function(rtab) rtab %>% melt(id.vars = "Q", value.name = "r") %>%
  mutate(variable = sub("^R_", "", variable))
rtab_indiv_melt <- melt_rtab(rtab_indiv)
rtab_rep_melt <- melt_rtab(rtab_rep)
rtab_rep_avg_melt <- melt_rtab(rtab_rep_avg)

# Make plots
plot_rtab <- function(rtab_melt, ylabel) ggplot(rtab_melt) + 
  geom_line(aes(x=Q, y=r, colour = variable), size = 1.5) +
  xlab("Diversity order (q)") + ylim(c(-1,0)) + ylab(ylabel) +
  scale_colour_discrete(name = "Metric",
                        labels = c("P20", "Zipf exponent (all)", 
                                   "Zipf exponent (filtered)")) +
  theme_classic() + theme_base

rplot_indiv <- plot_rtab(rtab_indiv_melt,
  expression("Pearson correlation with"~scriptstyle(""[q]*D^alpha)))
rplot_rep <- plot_rtab(rtab_rep_melt,
  expression("Pearson correlation with"~""[q]*D))
rplot_rep_avg <- plot_rtab(rtab_rep_avg_melt,
  expression("Pearson correlation with"~scriptstyle(""[q]*D^alpha)))
                      
rplot_out <- gplot_grid_onelegend(rplot_rep, rplot_rep_avg, nrow = 1,
                                  plot_height = plot_height, 
                                  plot_width = plot_width)

# Save optimum predictors
rbest_rep <- rtab_rep_melt %>% filter(r == min(r))
savetxt(rbest_rep$variable, outpath_rep_best_metric)
savetxt(rbest_rep$Q, outpath_rep_best_q)
savetxt(rbest_rep$r %>% round(2), outpath_rep_best_r)

rbest_avg <- rtab_rep_avg_melt %>% filter(r == min(r))
savetxt(rbest_avg$variable, outpath_avg_best_metric)
savetxt(rbest_avg$Q, outpath_avg_best_q)
savetxt(rbest_avg$r %>% round(2), outpath_avg_best_r)

rcross_avg <- rtab_rep_avg %>% filter(abs(R_P20) < abs(R_S_Filtered)) %>%
  filter(Q == min(Q))
savetxt(rcross_avg$Q, outpath_avg_cross_q)

rbest_avg_sfilter <- rtab_rep_avg_melt %>% filter(variable == "S_Filtered") %>%
  filter(r == min(r))
savetxt(rbest_avg_sfilter$Q, outpath_avg_sfilter_best_q)
savetxt(rbest_avg_sfilter$r %>% round(2), outpath_avg_sfilter_best_r)

#------------------------------------------------------------------------------
# SAVE FIGURES
#------------------------------------------------------------------------------
ggsave(plot = spectra, filename = outpath_spectra_grouped, device = "svg", units = "cm",
       height=plot_height, width = plot_width)
ggsave(plot = g_solo, filename = outpath_spectra_solo, device = "svg", units = "cm",
       height=15, width = 20)
ggsave(plot = g_indiv_alpha, filename = outpath_box_solo, device = "svg", units = "cm",
       height=20, width = 20*1.5)
ggsave(plot = rplot_out, filename = outpath_metrics_cor, device = "svg", units = "cm",
       height=plot_height, width = plot_width)
