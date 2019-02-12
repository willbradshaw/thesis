###############################################################################
## AUX FILE                                                                  ##
## Process Change-O read and counts tables                                   ##
###############################################################################

source("aux/packages.R")

#------------------------------------------------------------------------------
# COUNTS TABLES
#------------------------------------------------------------------------------

# Specify stage order 
changeo_count_stages <- c("presto_raw", "presto_filter", "presto_mask", 
                          "presto_correct", "presto_consensus", "presto_merged", 
                          "presto_collapsed", "presto_split", "changeo_makedb", 
                          "changeo_filter", "changeo_clonotyped", 
                          "changeo_germlined",
                          "changeo_split_functional")
changeo_count_names <- c("Raw reads", "Filtering", "Mask primers", "Cluster UMIs",
                         "Consensus reads", "Merge pairs", "Collapse identical sequences",
                         "Discard singletons", "Assign VDJs", "Filter by V-score",
                         "Clonotyping", "Germline assignment", "Functional sequences")

# Import counts table
import_counts <- function(path){
  col <- cols(NSEQS = "i", NREADS = "i", NREADS_RAW = "i", 
              CLUSTER_BARCODES = "d", CLUSTER_SETS = "d", NREADS_PC = "d",
              .default = col_character())
  tab <- import_tsv(path, col) %>% 
    arrange(match(STAGE, changeo_count_stages)) %>%
    mutate(STAGE = factor(STAGE, levels = changeo_count_stages))
  return(tab)
}

# Define plotting components
theme_countplot <- theme_classic() + theme_base + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = font),
        axis.title.x = element_blank())

# Set up plotting components
countplot_base <- function(colour = "REPLICATE", group = "REPLICATE", 
                           x = "STAGE") {
  ggplot(tab, mapping=aes_string(colour = colour, group = group, x = x)) +
  theme_countplot
}
countline_abs <- function(counts) geom_line(data = counts, aes(y = NREADS))
countline_rel <- function(counts) geom_line(data = counts, aes(y = NREADS_PC))
countscale_x <- function(stage_breaks = changeo_count_stages,
                         stage_names = changeo_count_names,
                         stages_include = changeo_count_stages){
  scale_x_discrete(breaks = factor(stage_breaks, levels = stage_breaks), 
                   labels = stage_names, limits = stages_include)
}
countscale_colour <- function(palette = "", name = "Replicate"){
  if (all(palette == "")){
    scale_colour_discrete(name = name)
  } else {
    scale_colour_manual(values = palette, name = name)
  }
}

# Set up complete plotting functions
countplot_abs <- function(counts, colour = "REPLICATE", group = "REPLICATE", 
                          x = "STAGE", stage_breaks = changeo_count_stages,
                          stage_names = changeo_count_names,
                          stages_include = changeo_count_stages,
                          cname = "Replicate", palette = ""){
  countplot_base(colour, group, x) + countline_abs(counts) + 
    scale_y_continuous(name = "# Reads", limits = c(0, max(counts$NREADS))) +
    countscale_x(stage_breaks, stage_names, stages_include) +
    countscale_colour(palette = palette, name = cname)
}
countplot_rel <- function(counts, colour = "REPLICATE", group = "REPLICATE", 
                          x = "STAGE", stage_breaks = changeo_count_stages,
                          stage_names = changeo_count_names,
                          stages_include = changeo_count_stages,
                          cname = "Replicate", palette = ""){
  countplot_base(colour, group, x) + countline_rel(counts) + 
    scale_y_continuous(name = "% Reads", limits = c(0, 1),
                       labels = function(y) round(y*100)) +
    countscale_x(stage_breaks, stage_names, stages_include) +
    countscale_colour(palette = palette, name = cname)
}

# Combine relative and absolute plots into a single figure with common legend
countplot_all <- function(counts, colour = "REPLICATE", group = "REPLICATE", 
                          x = "STAGE", stage_breaks = changeo_count_stages,
                          stage_names = changeo_count_names,
                          stages_include = changeo_count_stages,
                          cname = "Replicate", palette = ""){
  # Make separate absolute and relative read survival plots
  absplot <- countplot_abs(counts, colour, group, x, stage_breaks, stage_names,
                           stages_include, cname, palette)
  relplot <- countplot_rel(counts, colour, group, x, stage_breaks, stage_names,
                           stages_include, cname, palette)
  # Extract legend from absplot
  g <- ggplotGrob(absplot + theme(legend.position = "bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  # Combine plots without legend
  plt <- plot_grid(absplot + theme(legend.position = "none"),
                   relplot + theme(legend.position = "none"), 
                   ncol = 2, nrow = 1, labels="AUTO",
                   label_fontfamily = titlefont, label_fontface = "plain",
                   label_size = fontsize_base * fontscale_label)
  combined <- arrangeGrob(plt,
                          legend,
                          ncol = 1, nrow = 2,
                          heights = unit.c(unit(1, "npc") - lheight, lheight))
  # return gtable invisibly
  invisible(combined)
}


# 
# # Make plots for each experiment
# 
# # Make boxplots of final read survival for each experiment
# surv_all <- bind_rows(counts) %>% filter(STAGE == "changeo_split_functional")
# experiment_boxplots <- surv_all %>% 
#   filter(EXPERIMENT %in% c("pilot", "ageing", "gut")) %>% ggplot() + 
#   geom_boxplot(aes(x=EXPERIMENT, y=NREADS_PC, colour = EXPERIMENT)) + 
#   xlab("Experiment") + ylab("% Read Survival") + theme_classic() +
#   theme(legend.position = "none") + ylim(c(0,1))
# experiment_violins <- surv_all %>% 
#   filter(EXPERIMENT %in% c("pilot", "ageing", "gut")) %>% ggplot() + 
#   geom_violin(aes(x=EXPERIMENT, y=NREADS_PC, colour = EXPERIMENT)) + 
#   xlab("Experiment") + ylab("% Read Survival") + theme_classic() +
#   theme(legend.position = "none") + ylim(c(0,1))
# 
# # Repeat boxplotting averaged over individual
# surv_all_indiv <- surv_all %>% group_by(EXPERIMENT, INDIVIDUAL) %>%
#   summarise(NSEQS = sum(NSEQS), NREADS = sum(NREADS), NREADS_RAW = sum(NREADS_RAW),
#             NREADS_PC = NREADS/NREADS_RAW)
# experiment_boxplots_indiv <- surv_all_indiv %>% 
#   filter(EXPERIMENT %in% c("pilot", "ageing", "gut")) %>% ggplot() + 
#   geom_boxplot(aes(x=EXPERIMENT, y=NREADS_PC, colour = EXPERIMENT)) + 
#   xlab("Experiment") + ylab("% Read Survival") + theme_classic() +
#   theme(legend.position = "none") + ylim(c(0,1))
# experiment_violins_indiv <- surv_all_indiv %>% 
#   filter(EXPERIMENT %in% c("pilot", "ageing", "gut")) %>% ggplot() + 
#   geom_violin(aes(x=EXPERIMENT, y=NREADS_PC, colour = EXPERIMENT)) + 
#   xlab("Experiment") + ylab("% Read Survival") + theme_classic() +
#   theme(legend.position = "none") + ylim(c(0,1))
# 
# # Test for a difference in read survival between experiments
# surv_pilot <- surv_all %>% filter(EXPERIMENT == "pilot") %>% pull(NREADS_PC)
# surv_ageing <- surv_all %>% filter(EXPERIMENT == "ageing") %>% pull(NREADS_PC)
# surv_gut <- surv_all %>% filter(EXPERIMENT == "gut") %>% pull(NREADS_PC)
# 
# mw_p_g <- wilcox.test(surv_pilot, surv_gut)
# mw_p_g <- wilcox.test(surv_pilot, surv_ageing)
# 
# 

#------------------------------------------------------------------------------
# FULL CHANGE-O TABLES
#------------------------------------------------------------------------------

import_tab <- function(path){
  col <- cols(FUNCTIONAL = "l", IN_FRAME = "l", STOP = "l",
              MUTATED_INVARIANT = "l", INDELS = "l",
              V_SEQ_START = "i", V_SEQ_LENGTH = "i", V_GERM_START_VDJ = "i",
              V_GERM_LENGTH_VDJ = "i", V_GERM_START_IMGT = "i",
              V_GERM_LENGTH_IMGT = "i", NP1_LENGTH = "i",
              D_SEQ_START  = "i", D_SEQ_LENGTH = "i", D_GERM_START = "i",
              D_GERM_LENGTH = "i", NP2_LENGTH = "i",
              J_SEQ_START  = "i", J_SEQ_LENGTH = "i", J_GERM_START = "i",
              J_GERM_LENGTH = "i", JUNCTION_LENGTH = "i",
              V_SCORE = "d", V_IDENTITY = "d", V_EVALUE = "d", 
              D_SCORE = "d", D_IDENTITY = "d", D_EVALUE = "d", 
              J_SCORE = "d", J_IDENTITY = "d", J_EVALUE = "d",
              CONSCOUNT = "i", DUPCOUNT = "i",
              HAS_V = "l", HAS_D = "l", HAS_J = "l", 
              N_V = "i", N_D = "i", N_J = "i", V_AMBIG = "l", D_AMBIG = "l",
              J_AMBIG = "l", MU_FREQ_CDR_R = "d", MU_FREQ_CDR_S = "d",
              MU_FREQ_FWR_R = "d", MU_FREQ_FWR_S = "d",
              .default = col_character())
  tab <- import_tsv(path, col)
  return(tab)
}

diagnose_functionality <- function(tab){
  # Group reads by functionality and reason for nonfunctionality
  fail_descs <- list("000" = "Functional",
                     "100" = "No J",
                     "010" = "STOP",
                     "001" = "Frameshift",
                     "110" = "No J + STOP",
                     "101" = "No J + Frameshift",
                     "011" = "STOP + Frameshift",
                     "111" = "No J + STOP + Frameshift")
  tab <- tab %>%
    mutate(FUNC_STATE = paste0(as.integer(!HAS_J), as.integer(STOP), as.integer(!IN_FRAME)),
           FUNC_DESC = as.character(fail_descs[FUNC_STATE])
    )
  return(tab)
}


#------------------------------------------------------------------------------
# CLONE-SIZE DISTRIBUTIONS
#------------------------------------------------------------------------------

H <- function(n, s=1){
  # Compute the nth generalised harmonic number with exponent S
  sum(1/(1:n)^s)
}

lzipf <- function(f, s){
  # Compute the negative log-likelihood of a set of frequencies under a Zipf
  # distribution with parameter s
  s * sum(f * log(1:length(f))) + sum(f) * log(H(length(f), s))
}

dzipf <- function(r, s, N){
  # Return the predicted frequency of a rank r under a Zipf distribution
  # for a population with total size N
  (1/(r^s))/H(N, s)
}

compute_zipf_slope <- function(frequencies, n_exclude){
  m <- mle(function(s) lzipf(frequencies[-(1:n_exclude)], s), start = list(s=1))
  return(m@coef[["s"]])
}

compute_expected_frequencies_individual <- function(clntab, n_exclude){
  # Compute expected frequencies for a single individual with a single
  # exclusion number
  clntab %>% arrange(CLNRANK) %>% 
    mutate(S = compute_zipf_slope(CLNCOUNT, n_exclude),
           CLNFREQ_EXP = dzipf(CLNRANK, S, n()),
           CLNCOUNT_EXP = sum(CLNCOUNT)*CLNFREQ_EXP,
           CLNFREQ_RES = CLNFREQ/CLNFREQ_EXP,
           CHI_SQ = (CLNCOUNT - CLNCOUNT_EXP)^2/CLNCOUNT_EXP,
           N_EXCLUDE = n_exclude)
}

compute_expected_frequencies <- function(clntab, n_exclude){
  # Compute expected clonal frequencies under Zipf law using MLE
  # from a table of clonal ranks and frequencies
  tabs <- lapply(clntab %>% pull(INDIVIDUAL) %>% unique, function(i)
    compute_expected_frequencies_individual(filter(clntab, INDIVIDUAL == i),
                                            n_exclude[[i]]))
  return(bind_rows(tabs))
}

compute_clntab <- function(tab){
  tab %>% filter(!is.na(CLONE)) %>%
    group_by(INDIVIDUAL, CLONE) %>% summarise(CLNCOUNT = n()) %>%
    group_by(INDIVIDUAL) %>% arrange(CLNCOUNT) %>% 
    mutate(CLNRANK = row_number(desc(CLNCOUNT)),
           CLNFREQ = CLNCOUNT/sum(CLNCOUNT)) %>%
    arrange(CLNRANK)
}

plot_zipf_fit <- function(clntab){
  ggplot(clntab, aes(x=CLNRANK)) + 
    geom_point(aes(y=CLNFREQ, colour = INDIVIDUAL, 
                   shape = CLNRANK > N_EXCLUDE)) + 
    geom_line(aes(y=CLNFREQ_EXP)) + 
    scale_shape_manual(values = c(4, 16)) +
    scale_x_log10() + scale_y_log10() + 
    facet_wrap(~INDIVIDUAL, scales = "free") +
    ylab("Relative clonal frequency") +
    xlab("Clone rank in repertoire") +
    theme_classic() + theme_base +
    theme(legend.position = "none")
}

plot_zipf_residuals <- function(clntab){
  ggplot(clntab, aes(x=CLNRANK)) + 
    geom_point(aes(y=CLNFREQ/CLNFREQ_EXP, colour = INDIVIDUAL)) + 
    geom_hline(aes(yintercept = 1)) + 
    scale_x_log10() + scale_y_log10() + 
    facet_wrap(~INDIVIDUAL, scales = "free") +
    theme_classic() + theme_base
}

compute_chisq_p <- function(clntab, n_exclude){
  # Compute chi-squared goodness-of-fit p-values for the fit of a clone
  # table to a Zipf distribution
  tibble(INDIVIDUAL = names(n_exclude), N_EXCLUDE = n_exclude) %>% 
    inner_join(clntab, by="INDIVIDUAL") %>% filter(CLNRANK > N_EXCLUDE) %>%
    group_by(INDIVIDUAL) %>% summarise(CHI_SQ = sum(CHI_SQ), DF = n() - 1) %>% 
    mutate(P = pchisq(CHI_SQ, DF, lower.tail = FALSE))
}

#------------------------------------------------------------------------------
# DIVERSITY SPECTRA
#------------------------------------------------------------------------------

import_div <- function(path){
  col <- cols(Q = "d", N_GROUP = "i", D = "d", D_SD = "d", D_UPPER = "d",
              D_LOWER = "d", E = "d", E_UPPER = "d", E_LOWER = "d",
              .default = col_character())
  tab <- import_tsv(path, col)
  return(tab)
}

plot_diversity_spectra <- function(data, group_within, ...){
  ggplot(data) + 
    geom_line(aes_string(x="Q", y="D", colour = group_within)) + 
    geom_ribbon(aes_string(x="Q", ymin = "D_LOWER", ymax = "D_UPPER", 
                           fill = group_within, group = group_within), 
                alpha = 0.4) +
    xlab("Diversity order (q)") + 
    ylab(expression(Diversity~(""[q]*D))) +
    theme_classic() + theme_base
} # TODO: Add theme, axis titles etc.

plot_evenness_spectra <- function(data, group_within, ...){
  ggplot(data) + 
    geom_line(aes_string(x="Q", y="E", colour = group_within)) + 
    geom_ribbon(aes_string(x="Q", ymin = "E_LOWER", ymax = "E_UPPER", 
                           fill = group_within, group = group_within), 
                alpha = 0.4)
}

plot_diversity_gamma <- function(data, group_within){
  data %>% filter(DIVTYPE == "gamma") %>% plot_diversity_spectra(group_within) +
    ylab(expression("Gamma-diversity"~(scriptstyle(""[q]*D^gamma))))
}

plot_diversity_alpha <- function(data, group_within){
  data %>% filter(DIVTYPE == "alpha") %>% plot_diversity_spectra(group_within) +
    ylab(expression("Alpha-diversity"~(scriptstyle(""[q]*D^alpha))))
}

plot_diversity_beta <- function(data, group_within){
  data %>% filter(DIVTYPE == "beta") %>% plot_diversity_spectra(group_within) +
    ylab(expression("Beta-diversity"~(scriptstyle(""[q]*D^beta))))
}


rescale_summarised_beta <- function(beta_div_tab){
  # Rescale beta-diversity to scale from 0 (min) to 1 (max)
  beta_div_tab %>%
    mutate(D = (D-1)/(N_GROUP-1),
           D_UPPER = (D_UPPER-1)/(N_GROUP-1),
           D_LOWER = (D_LOWER-1)/(N_GROUP-1),
           D_SD = D_SD/(N_GROUP-1)) %>%
    select(-N_GROUP)
}

plot_diversity_beta_scaled <- function(data, group_within){
  data %>% filter(DIVTYPE == "beta") %>% rescale_summarised_beta %>%
    plot_diversity_spectra(group_within) + ylim(c(NA,1)) +
    ylab(expression("Rescaled beta-diversity"~(scriptstyle(""[q]*D^beta))))
}

#------------------------------------------------------------------------------
# MODELLING SPECIFIC DIVERSITIES
#------------------------------------------------------------------------------




