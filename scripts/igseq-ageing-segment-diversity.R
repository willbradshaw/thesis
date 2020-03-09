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
# IMPORT DATA
#------------------------------------------------------------------------------

tab_grouped <- import_div(inpath_grouped)
tab_solo <- import_div(inpath_solo)

#------------------------------------------------------------------------------
# GENERATE ALPHA SPECTRA
#------------------------------------------------------------------------------

# Parameters
palette <- c(colours_igseq[["ageing_group1"]], colours_igseq[["ageing_group2"]],
             colours_igseq[["ageing_group3"]], colours_igseq[["ageing_group4"]])
# palette <- gg_color_hue(4)
# TODO: get age-group palette
age_groups <- c("39", "56", "73", "128")
tab_grouped <- tab_grouped %>% mutate(AGE_DAYS = factor(AGE_DAYS, levels = age_groups))
tab_solo <- tab_solo %>% mutate(AGE_DAYS = factor(AGE_DAYS, levels = age_groups))

g_alpha <- plot_diversity_alpha(tab_grouped, "AGE_DAYS") +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm")),
        legend.text = element_text(margin=margin(r = 0.3, unit="cm")))
g_beta_unscaled <- plot_diversity_beta(tab_grouped, "AGE_DAYS") +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm")),
        legend.text = element_text(margin=margin(r = 0.3, unit="cm")))
g_beta_scaled <- plot_diversity_beta_scaled(tab_grouped, "AGE_DAYS") +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  ylim(c(0,0.25)) +
  theme(legend.title = element_text(margin = margin(r = 1, unit = "cm")),
        legend.text = element_text(margin=margin(r = 0.3, unit="cm")))


g_solo <-  ggplot(tab_solo) + 
  geom_line(aes(x=Q, y=D, colour = AGE_DAYS, group = INDIVIDUAL)) + 
  geom_ribbon(aes(x=Q, ymin = D_LOWER, ymax = D_UPPER, 
                  fill = AGE_DAYS, group = INDIVIDUAL), alpha = 0.4) +
  facet_wrap(~AGE_DAYS, scales = "free") +
  xlab("Diversity order (q)") + 
  ylab(expression(Diversity~(""[q]*D))) +
  xlim(c(0,4)) + ylim(c(0,120)) +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  theme_classic() + theme_base

#------------------------------------------------------------------------------
# COMBINE ALPHA AND BETA SPECTRA WITH SINGLE LEGEND
#------------------------------------------------------------------------------

# Extract legend from absplot
g <- ggplotGrob(g_alpha)$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$width)
# Combine plots without legend
plt <- plot_grid(g_alpha + theme(legend.position = "none"),
                 g_beta_scaled + theme(legend.position = "none"), 
                 ncol = 2, nrow = 1, labels="AUTO",
                 label_fontfamily = titlefont, label_fontface = "plain",
                 label_size = fontsize_base * fontscale_label)
combined <- arrangeGrob(plt,
                        legend,
                        ncol = 1, nrow = 2,
                        heights = unit.c(unit(1, "npc") - lheight, lheight))

# Visualise plot
plot_unit = "cm"
plot_height <- 15
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

plt <- grid.grab()

#------------------------------------------------------------------------------
# SAVE SPECTRA
#------------------------------------------------------------------------------

ggsave(plot = plt, filename = outpath_alpha_beta, device = "svg", units = "cm",
       height=15, width = 25)

#------------------------------------------------------------------------------
# COMPARE DIFFERENT AGE GROUPS AT SPECIFIC ORDERS
#------------------------------------------------------------------------------

# Filter solo-diversity table
filter_solo_tab <- function(q){
  tab_solo %>% filter(Q == q) %>%
    select(AGE_DAYS, INDIVIDUAL, D) %>%
    mutate(AGE_DAYS = as.numeric(levels(AGE_DAYS))[AGE_DAYS]) %>%
    arrange(AGE_DAYS, INDIVIDUAL)
}
multi_filter <- function(qvals, family = Gamma()){
  bind_rows(lapply(qvals, function(q) 
    filter_solo_tab(q) %>% mutate(Q = q)))
}

# Generate (G)LMs
glm_diversity <- function(q, family = Gamma()){
  if (all(is.na(family))){
    lm(D~AGE_DAYS, data = filter_solo_tab(q))
  } else {
    glm(D~AGE_DAYS, family = family,
        data = filter_solo_tab(q))
  }
}
predict_glm_diversity <- function(q, family = Gamma()){
  g <- glm_diversity(q, family)
  ages <- tibble(AGE_DAYS = seq(min(as.numeric(age_groups)),
                                max(as.numeric(age_groups))))
  if (all(is.na(family))){
    p <- predict.lm(g, ages, type = "response")
  } else {
    p <- predict.glm(g, ages, type = "response")
  }
  return(bind_cols(ages, D = p))
}

multi_predict <- function(qvals, family = Gamma()){
  bind_rows(lapply(qvals, function(q) 
    predict_glm_diversity(q, family) %>% mutate(Q = q)))
}

# Kruskal-Wallis non-parametric ANOVA
kruskal_diversity <- function(q){
  kruskal.test(formula = D ~ AGE_DAYS, data = filter_solo_tab(q))
}

# Generate summary tables
summarise_glm_diversity <- function(q, family = Gamma()){
  l <- glm_diversity(q, family)
  s <- summary(l)
  return(s)
}
signif_table <- function(qvals, families = list("gamma" = Gamma())){
  # Compute p-values under GLM families
  L <- lapply(families, function(f) sapply(qvals, function(q)
    summary(glm_diversity(q,f))$coefficients[2,4]))
  # Compute p-values under Kruskal-Wallis ANOVA
  K <- lapply(qvals, function(q) kruskal_diversity(q)$p.value)
  # Combine together
  S_GLM <- melt(L) %>% rename(FAMILY = L1, P = value) %>% 
    group_by(FAMILY) %>% mutate(Q = qvals)
  S_Kruskal <- melt(K) %>% rename(P = value) %>%
    mutate(FAMILY = "kruskal", Q = qvals[L1]) %>% select(-L1)
  return(bind_rows(S_GLM, S_Kruskal))
}

# Generate plots
plot_solo_diversity <- function(qvals, family = Gamma()){
  ggplot(multi_filter(qvals), aes(x=AGE_DAYS, y=D)) +
    geom_boxplot(aes(fill = factor(AGE_DAYS, levels = age_groups),
                     group = AGE_DAYS), outlier.shape = NA) +
    geom_point(alpha = 0.4, size = 2) +
    geom_line(data = multi_predict(qvals, family), colour = "black", 
              size = 2) +
    facet_wrap(~Q, scales = "free", labeller = function(q) label_both(q, sep = " = ")) +
    scale_fill_manual(values = palette, name = "Age group (days)") +
    ylab("Diversity") + xlab("Age at death (days)") +
    theme_classic() + theme_base
}
annotate_pvalues <- function(qvals, family = list("gamma" = Gamma()),
                             nsig = 2){
  S <- signif_table(qvals, family) %>% mutate(P = signif(P, nsig))
  S_GLM <- filter(S, FAMILY == names(family))
  S_KWT <- filter(S, FAMILY == "kruskal")
  T <- tibble(Q = qvals) %>% 
    mutate(ANNOT_GLM = paste0("P(GLM) = ", filter(S_GLM, Q == Q) %>% pull(P)),
           ANNOT_KWT = paste0("P(KWT) = ", filter(S_KWT, Q == Q) %>% pull(P)),
           ANNOT = paste(ANNOT_GLM, ANNOT_KWT, sep = "\n"))
  return(T)
}

# Make boxplots and annotate with P-values
qvals <- c(0,1,1.5,2,3,4)
P_pos <- tibble(Q = qvals, x = 110, y = c(97.5,42,28.5,21,14,11.5))

P_gamma <- annotate_pvalues(qvals) %>% full_join(P_pos, by = "Q")
P_linear <- annotate_pvalues(qvals, list("linear" = NA)) %>%
  full_join(P_pos, by = "Q")
P_igauss <- annotate_pvalues(qvals, list("igauss" = inverse.gaussian())) %>%
  full_join(P_pos, by = "Q")

g_solofit_gamma <- plot_solo_diversity(qvals, family = Gamma()) +
  geom_text(data = P_gamma, aes(label = ANNOT,x=x,y=y), size = 3.5)
g_solofit_linear <- plot_solo_diversity(qvals, family = NA) +
  geom_text(data = P_linear, aes(label = ANNOT,x=x,y=y), size = 3.5)
g_solofit_igauss <- plot_solo_diversity(qvals, family = inverse.gaussian()) +
  geom_text(data = P_igauss, aes(label = ANNOT,x=x,y=y), size = 3.5)


#------------------------------------------------------------------------------
# SAVE FIGURES
#------------------------------------------------------------------------------

ggsave(plot = g_alpha, filename = outpath_spectra_grouped, device = "svg", units = "cm",
       height=plot_height, width = plot_width)
ggsave(plot = g_solo, filename = outpath_spectra_solo, device = "svg", units = "cm",
       height=15, width = 20)
ggsave(plot = g_solofit_gamma, filename = outpath_solofit_gamma, device = "svg", units = "cm",
       height=20, width = 20*1.5)
ggsave(plot = g_solofit_linear, filename = outpath_solofit_linear, device = "svg", units = "cm",
       height=20, width = 20*1.5)
ggsave(plot = g_solofit_igauss, filename = outpath_solofit_igauss, device = "svg", units = "cm",
       height=20, width = 20*1.5)
