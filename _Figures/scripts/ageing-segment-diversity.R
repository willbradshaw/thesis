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
segments <- "VJ" # or "VJ"
tab_path_grouped <- paste0("../_Data/changeo/spectra/ageing_", segments,
                           "-diversity-grouped_seqs-",
                           seqset, "_copy-", copy, ".tsv")
tab_path_solo <- paste0("../_Data/changeo/spectra/ageing_", segments,
                        "-diversity-solo_seqs-",
                        seqset, "_copy-", copy, ".tsv")


# Output paths
filename_base <- paste0("ageing-", segments, "-diversity")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab_grouped <- import_div(tab_path_grouped)
tab_solo <- import_div(tab_path_solo)

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
  scale_fill_manual(values = palette, name = "Age group (days)")
g_beta_unscaled <- plot_diversity_beta(tab_grouped, "AGE_DAYS") +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_fill_manual(values = palette, name = "Age group (days)")
g_beta_scaled <- plot_diversity_beta_scaled(tab_grouped, "AGE_DAYS") +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  ylim(c(0,0.3))

#------------------------------------------------------------------------------
# COMPARE SHANNON ENTROPY OF DIFFERENT AGE GROUPS
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
                     group = AGE_DAYS)) +
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

savefig(plot = g_alpha, filename = paste0(filename_base, "-alpha"),
        height = 20, ratio = 1.5)
savefig(plot = g_beta_unscaled, filename = paste0(filename_base, "-beta-unscaled"),
        height = 20, ratio = 1.5)
savefig(plot = g_beta_scaled, filename = paste0(filename_base, "-beta-scaled"),
        height = 20, ratio = 1.5)
savefig(plot = g_solofit_gamma, height = 20, ratio = 1.5,
        filename = paste0(filename_base, "-solo-fit-gamma"))
savefig(plot = g_solofit_linear, height = 20, ratio = 1.5,
        filename = paste0(filename_base, "-solo-fit-linear"))
savefig(plot = g_solofit_igauss, height = 20, ratio = 1.5,
        filename = paste0(filename_base, "-solo-fit-igauss"))