###############################################################################
## FIGURE                                                                    ##
## CLONE DISTRIBUTION CORRELATION BETWEEN REPLICATES                         ##
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

# Input paths
tab_path <- "../_Data/changeo/ctabs/pilot-final.tab"

# Output paths
filename_nclones <- "pilot-nclones"
filename_base <- "pilot-clone-sizes"

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(tab_path)

#------------------------------------------------------------------------------
# COUNT UNIQUE SEQUENCES WITH NA CLONES
#------------------------------------------------------------------------------

nseq_naclone <- tab %>% group_by(SEQUENCE_INPUT, CLONE) %>% summarise() %>% 
  pull(CLONE) %>% is.na %>% mean %>% (function(x) round((1-x)*100, 1))
savetxt(nseq_naclone, "pilot-nseq-assigned-clones")

#------------------------------------------------------------------------------
# COUNT CLONES PER INDIVIDUAL
#------------------------------------------------------------------------------

tab_cl <- tab %>% mutate(REP = sub("\\d-\\d\\d", "", REPLICATE)) %>%
  filter(!is.na(CLONE)) %>%
  group_by(INDIVIDUAL, REP, CLONE) %>% 
  summarise(CLNCOUNT = n(), DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT))

tab_cl_counts <- tab_cl %>% group_by(INDIVIDUAL, CLONE) %>% 
  summarise(CLNCOUNT = sum(CLNCOUNT), DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT)) %>%
  group_by(INDIVIDUAL) %>% 
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

#------------------------------------------------------------------------------
# COUNT CLONES BY SIZE
#------------------------------------------------------------------------------

tab_cl_size <- tab_cl %>% group_by(INDIVIDUAL, CLONE) %>% 
  summarise(CLNCOUNT = sum(CLNCOUNT)) %>% group_by(CLNCOUNT) %>% 
  summarise(N = n()) %>% mutate(N_PC = N/sum(N) * 100,
                                CUM_PC = cumsum(N_PC))

nclones_1count <- tab_cl_size %>% filter(CLNCOUNT == 1) %>% pull(N_PC) %>%
  round(1)
nclones_small <- tab_cl_size %>% filter(CLNCOUNT < 5) %>% pull(N_PC) %>%
  sum %>% round(1)
nclones_large <- tab_cl_size %>% filter(CLNCOUNT >= 5) %>% pull(N_PC) %>%
  sum %>% round(1)
savetxt(nclones_1count, paste0(filename_nclones, "-pc-1count"))
savetxt(nclones_small, paste0(filename_nclones, "-pc-small"))

#------------------------------------------------------------------------------
# COUNT CLONES PRESENT IN 1/2/3 replicates
#------------------------------------------------------------------------------

tab_cl_spread <- spread(tab_cl, REP, CLNCOUNT) %>%
  mutate(bio = ifelse(is.na(bio), 0, bio),
         lib = ifelse(is.na(lib), 0, lib),
         orig = ifelse(is.na(orig), 0, orig)) %>%
  group_by(INDIVIDUAL, CLONE) %>%
  summarise(DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT),
            bio = sum(bio), lib = sum(lib), orig = sum(orig)) %>%
  mutate(N_ABSENT = ((bio==0) + (lib==0) + (orig == 0)),
         N_PRESENT = 3-N_ABSENT,
         CLNCOUNT = (bio+lib+orig),
         CLNCOUNT_AVG = (bio+lib+orig)/(N_PRESENT))

tab_cl_nrep <- tab_cl_spread %>% 
  mutate(N_ABSENT = ((bio==0) + (lib==0) + (orig == 0)),
         N_PRESENT = 3-N_ABSENT)

tab_cl_nrep_summ <- tab_cl_nrep %>% group_by(N_PRESENT) %>% 
  summarise(N = n()) %>% mutate(N_PC = N/sum(N)*100)

# Extract text values
nclones_1rep <- tab_cl_nrep_summ %>% filter(N_PRESENT == 1) %>%
  pull(N_PC) %>% round(2)
nclones_2rep <- tab_cl_nrep_summ %>% filter(N_PRESENT == 2) %>%
  pull(N_PC) %>% round(2)
nclones_3rep <- tab_cl_nrep_summ %>% filter(N_PRESENT == 3) %>%
  pull(N_PC) %>% round(2)

savetxt(nclones_1rep, paste0(filename_nclones, "-pc-1rep"))
savetxt(nclones_2rep, paste0(filename_nclones, "-pc-2rep"))
savetxt(nclones_3rep, paste0(filename_nclones, "-pc-3rep"))

#------------------------------------------------------------------------------
# PLOT CLONE-SIZE DISTRIBUTION
#------------------------------------------------------------------------------

palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
             colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]])

g_size <- tab_cl %>% group_by(INDIVIDUAL, CLONE) %>% 
  summarise(CLNCOUNT = sum(CLNCOUNT)) %>% group_by(INDIVIDUAL, CLNCOUNT) %>% 
  summarise(N = n()) %>% mutate(N_PC = N/sum(N) * 100,
                                CUM_PC = cumsum(N_PC)) %>%
  ggplot() + geom_line(aes(x=CLNCOUNT, y=N_PC, colour = INDIVIDUAL), 
                       size = 1.5) +
  xlab("# Unique sequences per clone") + ylab("% of clones") +
  scale_colour_manual(values = palette, name = "Individual") +
  scale_x_log10() + theme_classic() + theme_base +
  theme(legend.title = element_text(margin = margin(r=0.5, unit = "cm")))

#------------------------------------------------------------------------------
# PLOT CLNCOUNT/NREP RELATIONSHIP
#------------------------------------------------------------------------------

g_nrep <- tab_cl_nrep %>% group_by(CLNCOUNT, N_PRESENT) %>% 
  summarise(N = n()) %>% group_by(CLNCOUNT) %>% mutate(N_PC = N/sum(N)) %>%
  ggplot() + geom_smooth(aes(x=CLNCOUNT, y=N_PC,
                              colour = as.factor(N_PRESENT),
                              fill = as.factor(N_PRESENT)), 
                          method = "loess", formula = y~x) + 
  scale_x_log10(name = "# Unique sequences per clone") + 
  scale_y_continuous(name = "% of clones", labels = function(x) round(x*100)) +
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

#------------------------------------------------------------------------------
# GET CORRELATION IN CLONE SIZES
#------------------------------------------------------------------------------

tab_cl_cor <- tab_cl_spread %>%
  summarise(BL = cor(bio, lib), BO = cor(bio, orig), LO = cor(lib, orig)) %>%
  melt(id.vars = "INDIVIDUAL", variable.name = "COMPARISON", 
       value.name = "R")

interrep_cor_avg <- tab_cl_cor %>% pull(R) %>% mean %>% round(2)
savetxt(interrep_cor_avg, paste0(filename_base, "-cor-avg"))

# Plot distribution of inter-replicate correlation coefficients
g_cor <- ggplot(tab_cl_cor, aes(x=COMPARISON, y=R)) + 
  geom_boxplot(aes(fill = COMPARISON), outlier.shape = NA) +
  geom_point(size = 3, alpha = 0.4) +
  scale_x_discrete(labels = c("bio/lib", "bio/orig", "lib/orig"),
                   name = "Inter-replicate comparison") + 
  scale_y_continuous(breaks = seq(0,1,0.2),
                   limits = c(0,1), name = "Correlation coefficient (r)") + 
  theme_classic() + theme_base +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = fontsize_base))

savefig(g_cor, paste0(filename_base, "-cor-boxplots"), 
        height = 10*1.2, width = 15*1.2)

#------------------------------------------------------------------------------
# PLOT INTER-REPLICATE CORRELATIONS ON SCATTER PLOTS
#------------------------------------------------------------------------------

g_inter <- ggplot(tab_cl_spread) + 
  geom_point(aes(x=orig, y=bio, colour = INDIVIDUAL)) +
  facet_wrap(~INDIVIDUAL, scales = "free") +
  scale_colour_manual(values = palette, name = "Individual") +
  xlab("Replicate 1 (orig)") + ylab("Replicate 2 (bio)") +
  scale_x_log10() + scale_y_log10() +
  theme_classic() + theme_base + theme(legend.position = "none")

savefig(g_inter, paste0(filename_base, "-cor-scatter"), 
        height = 20, width = 20)

  