###############################################################################
## FIGURE                                                                    ##
## Clone counts and statistics in ageing dataset                             ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

#rm(list = ls())

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")
source("aux/ggplot2.R")
source("aux/changeo.R")

# Input paths
tab_path <- "../_Data/changeo/ctabs/ageing-final.tab"

# Output paths
filename_nclones <- "ageing-nclones"
filename_base <- "ageing-clone-sizes"

# Parameters
palette <- c(colours_igseq[["ageing_group1"]], colours_igseq[["ageing_group2"]],
             colours_igseq[["ageing_group3"]], colours_igseq[["ageing_group4"]])
age_groups <- c("39", "56", "73", "128")

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(tab_path)

#------------------------------------------------------------------------------
# COUNT UNIQUE SEQUENCES WITH NA CLONES
#------------------------------------------------------------------------------

nseq_naclone <- tab %>% group_by(SEQUENCE_INPUT, CLONE) %>% summarise() %>% 
  pull(CLONE) %>% is.na %>% mean %>% (function(x) round((1-x)*100, 1))
savetxt(nseq_naclone, "ageing-nseq-assigned-clones")

#------------------------------------------------------------------------------
# COUNT CLONES PER INDIVIDUAL
#------------------------------------------------------------------------------

tab_cl <- tab %>% mutate(REP = sub("\\d-\\d\\d", "", REPLICATE)) %>%
  filter(!is.na(CLONE)) %>%
  group_by(AGE_DAYS, INDIVIDUAL, REP, CLONE) %>% 
  summarise(CLNCOUNT = n(), DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT))

tab_cl_counts <- tab_cl %>% group_by(AGE_DAYS, INDIVIDUAL, CLONE) %>% 
  summarise(CLNCOUNT = sum(CLNCOUNT), DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT)) %>%
  group_by(AGE_DAYS, INDIVIDUAL) %>% 
  summarise(N = n(), CLNCOUNT = sum(CLNCOUNT), DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT))

# Extract and save text values
clones_individual_min <- tab_cl_counts %>% pull(N) %>% min %>% 
  (function(x) floor(x/10^2)*10^2)
clones_individual_max <- tab_cl_counts %>% pull(N) %>% max %>% 
  (function(x) ceiling(x/10^2)*10^2)
clones_individual_med <- tab_cl_counts %>% pull(N) %>% median %>% round
savetxt(clones_individual_min, paste0(filename_nclones, "-individual-min"))
savetxt(clones_individual_max, paste0(filename_nclones, "-individual-max"))
savetxt(clones_individual_med, paste0(filename_nclones, "-individual-med"))

# Counts of individuals from pilot study
tab_cl_counts_pilot_indivs <- tab_cl_counts %>% 
  filter(INDIVIDUAL %in% c("2-03", "2-04", "2-05", "2-06"))
clones_individual_min_pilot <- tab_cl_counts_pilot_indivs %>% pull(N) %>% min %>% 
  (function(x) floor(x/10^2)*10^2)
clones_individual_max_pilot <- tab_cl_counts_pilot_indivs %>% pull(N) %>% max %>% 
  (function(x) ceiling(x/10^2)*10^2)
savetxt(clones_individual_min_pilot, paste0(filename_nclones, "-individual-min-pilot"))
savetxt(clones_individual_max_pilot, paste0(filename_nclones, "-individual-max-pilot"))

# Test age effect on clone counts
tab_cl_counts_num <- tab_cl_counts %>% ungroup() %>% 
  mutate(AGE_DAYS = as.numeric(AGE_DAYS))
cl_counts_kruskal <- kruskal.test(formula = N ~ AGE_DAYS, data = tab_cl_counts_num)
savetxt(cl_counts_kruskal$p.value %>% signif(2), 
        paste0(filename_nclones, "-kruskal-p"))

# Plot range of clone numbers in different age groups
age_groups <- c("39", "56", "73", "128")
g_nclones <- ggplot(tab_cl_counts_num) + 
  geom_boxplot(aes(x=AGE_DAYS, y=N, fill = factor(AGE_DAYS, levels = age_groups))) +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  xlab("Age at death (days)") + ylab("# Clones") +
  theme_classic() + theme_base
savefig(g_nclones, filename_nclones, height = 15, ratio = 1.5)

#------------------------------------------------------------------------------
# COUNT CLONES BY SIZE
#------------------------------------------------------------------------------

tab_cl_size <- tab_cl %>% group_by(AGE_DAYS, INDIVIDUAL, CLONE) %>% 
  summarise(CLNCOUNT = sum(CLNCOUNT), DUPCOUNT = sum(DUPCOUNT),
            CONSCOUNT = sum(CONSCOUNT)) %>% group_by(CLNCOUNT) %>% 
  summarise(N = n(), DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT)) %>% 
  mutate(N_PC = N/sum(N) * 100, CUM_PC = cumsum(N_PC),
         CLNCOUNT_TOTAL = CLNCOUNT * N, 
         CLNCOUNT_PC = CLNCOUNT_TOTAL/sum(CLNCOUNT_TOTAL) * 100)

# Number of small clones
nclones_1count <- tab_cl_size %>% filter(CLNCOUNT == 1) %>% pull(N_PC) %>%
  round(1)
nclones_small <- tab_cl_size %>% filter(CLNCOUNT < 5) %>% pull(N_PC) %>%
  sum %>% round(1)
nclones_large <- tab_cl_size %>% filter(CLNCOUNT >= 5) %>% pull(N_PC) %>%
  sum %>% round(1)
savetxt(nclones_1count, paste0(filename_nclones, "-pc-1count"))
savetxt(nclones_small, paste0(filename_nclones, "-pc-small"))

#------------------------------------------------------------------------------
# COUNT NUMBER OF SEQUENCES IN EACH SIZE CLONE IN EACH AGE GROUP
#------------------------------------------------------------------------------

# Get clone-size counts for each individual
tab_cl_size_age <- tab_cl %>% group_by(AGE_DAYS, INDIVIDUAL, CLONE) %>% 
  summarise(CLNCOUNT = sum(CLNCOUNT), DUPCOUNT = sum(DUPCOUNT),
            CONSCOUNT = sum(CONSCOUNT)) %>% group_by(AGE_DAYS, INDIVIDUAL, CLNCOUNT) %>% 
  summarise(N = n(), DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT)) %>% 
  mutate(N_PC = N/sum(N) * 100, CUM_PC = cumsum(N_PC),
         CLNCOUNT_TOTAL = CLNCOUNT * N, 
         CLNCOUNT_PC = CLNCOUNT_TOTAL/sum(CLNCOUNT_TOTAL) * 100)

# Compute average proportion of abundant vs non-abundant by age group
tab_cl_abundant_indiv <- tab_cl_size_age %>% 
  group_by(AGE_DAYS, INDIVIDUAL, ABUNDANT = CLNCOUNT >= 5) %>%
  summarise(N_PC = sum(N_PC), CLNCOUNT_PC = sum(CLNCOUNT_PC))
tab_cl_abundant_age <- tab_cl_abundant_indiv %>%
  group_by(AGE_DAYS, ABUNDANT) %>% 
  summarise(N_PC = mean(N_PC), CLNCOUNT_PC = mean(CLNCOUNT_PC))

# Extract and save proportion of repertoire sequences contained in small clones
nseq_small_avg <- tab_cl_abundant_indiv %>% filter(!ABUNDANT) %>% 
  pull(CLNCOUNT_PC) %>% mean %>% round(1)
savetxt(nseq_small_avg, "ageing-pc-seq-in-small-clones-avg")

# Test age effect on abundance
nseq_age_kruskal <- tab_cl_abundant_indiv %>% filter(!ABUNDANT) %>%
  kruskal.test(formula = CLNCOUNT_PC ~ as.numeric(AGE_DAYS), data = .)
savetxt(nseq_age_kruskal$p.value %>% signif(2), 
        paste0("ageing-pc-seq-in-small-clones-kruskal-p"))

# Plot change in abundance proportion with age
g_nseq_abundant <- tab_cl_abundant_age %>% ggplot() +
  geom_col(aes(x=factor(AGE_DAYS, levels = age_groups), 
               y=CLNCOUNT_PC, fill=factor(AGE_DAYS, levels = age_groups), 
               alpha=factor(ABUNDANT, levels = c(TRUE, FALSE)))) +
  scale_fill_manual(values = palette, name = "Age group (days)") +
  scale_alpha_manual(values = c(0.3, 1), name = "Sequences in clone",
                     labels = c("5 or more", "4 or fewer")) +
  xlab("Age at death (days)") + ylab("Proportion of unique sequences (%)") +
  theme_classic() + theme_base + theme(
    legend.box = "vertical",
    legend.margin = margin(b=0.05, unit = "cm"),
    legend.box.margin = margin(t=0.2, unit = "cm")
  )
savefig(g_nseq_abundant, "ageing-pc-seq-in-small-clones", height = 15, ratio = 1.5)

#------------------------------------------------------------------------------
# COUNT CLONES PRESENT IN 1/2/3 replicates
#------------------------------------------------------------------------------

tab_cl_spread <- tab_cl %>% group_by(AGE_DAYS, INDIVIDUAL) %>%
  spread(REP, CLNCOUNT) %>% group_by(AGE_DAYS, INDIVIDUAL, CLONE) %>%
  mutate(A = ifelse(is.na(A), 0, A),
         B = ifelse(is.na(B), 0, B)) %>%
  summarise(DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT),
            A = sum(A), B = sum(B)) %>%
  mutate(N_ABSENT = ((A==0) + (B==0)),
         N_PRESENT = 2-N_ABSENT,
         CLNCOUNT = (A+B),
         CLNCOUNT_AVG = (A+B)/(N_PRESENT))

tab_cl_nrep_summ <- tab_cl_spread %>% group_by(N_PRESENT) %>% 
  summarise(N = n()) %>% mutate(N_PC = N/sum(N)*100)

# Extract text values
nclones_1rep <- tab_cl_nrep_summ %>% filter(N_PRESENT == 1) %>%
  pull(N_PC) %>% round(2)
nclones_2rep <- tab_cl_nrep_summ %>% filter(N_PRESENT == 2) %>%
  pull(N_PC) %>% round(2)

savetxt(nclones_1rep, paste0(filename_nclones, "-pc-1rep"))
savetxt(nclones_2rep, paste0(filename_nclones, "-pc-2rep"))

#------------------------------------------------------------------------------
# PLOT CLONE-SIZE DISTRIBUTION
#------------------------------------------------------------------------------

g_size <- tab_cl %>% group_by(AGE_DAYS, CLONE) %>% 
  summarise(CLNCOUNT = sum(CLNCOUNT)) %>% group_by(AGE_DAYS, CLNCOUNT) %>% 
  summarise(N = n()) %>% mutate(N_PC = N/sum(N) * 100,
                                CUM_PC = cumsum(N_PC)) %>%
  ggplot() + geom_line(aes(x=CLNCOUNT, y=N_PC, 
                           colour = factor(AGE_DAYS, levels = age_groups)), 
                       size = 1.5) +
  xlab("# Unique sequences per clone") + ylab("% of clones") +
  scale_colour_manual(values = palette, name = "Age group (days)") +
  scale_x_log10() + theme_classic() + theme_base +
  theme(legend.title = element_text(margin = margin(r=1, unit = "cm")))

#------------------------------------------------------------------------------
# PLOT CLNCOUNT/NREP RELATIONSHIP
#------------------------------------------------------------------------------

g_nrep <- tab_cl_spread %>% group_by(CLNCOUNT, N_PRESENT) %>% 
  summarise(N = n()) %>% group_by(CLNCOUNT) %>% mutate(N_PC = N/sum(N)) %>%
  ggplot() + geom_smooth(aes(x=CLNCOUNT, y=N_PC,
                              colour = as.factor(N_PRESENT),
                              fill = as.factor(N_PRESENT)), 
                          formula = y~x, method = "loess") + 
  scale_x_log10(name = "# Unique sequences per clone") + ylab("% of clones") +
  scale_fill_discrete(name = "Number of Replicates") + 
  scale_colour_discrete(name = "Number of Replicates") +
  coord_cartesian(ylim=c(0,1)) + theme_classic() + theme_base +
  theme(legend.title = element_text(margin = margin(r=1, unit = "cm")))


#------------------------------------------------------------------------------
# ARRANGE AND SAVE COUNT PLOTS
#------------------------------------------------------------------------------

# TODO: Increase spacing between legend title and labels
plt <- plot_grid(g_size, g_nrep,
          ncol = 2, nrow = 1, labels="AUTO",
          label_fontfamily = titlefont, label_fontface = "plain",
          label_size = fontsize_base * fontscale_label)

plot_height<- 10
plot_width <- 25
savefig(plt, filename_base, height = plot_height, width = plot_width)