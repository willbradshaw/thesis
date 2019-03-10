###############################################################################
## FIGURE                                                                    ##
## FITTING RANK/ABUNDANCE DISTRIBUTIONS OF AGEING CLONES TO ZIPF'S LAW       ##
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
tab_path <- "../_Data/changeo/ctabs/ageing-final.tab"

# Output paths
filename_base <- "ageing-clones-zipf"

# Set parameters
palette <- c(colours_igseq[["ageing_group1"]], colours_igseq[["ageing_group2"]],
             colours_igseq[["ageing_group3"]], colours_igseq[["ageing_group4"]])
n_exclude_value <- 5
zexp_round_ndigits <- 3
zexp_annot_x <- 800
zexp_annot_y <- 10^(-1.5)

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(tab_path)
age_groups <- tab %>% pull(AGE_DAYS) %>% unique %>% as.numeric %>% sort

#------------------------------------------------------------------------------
# COMPUTE CLONAL RANKS AND SIZES
#------------------------------------------------------------------------------

# Basic clone table
clntab_base <- tab %>% group_by(AGE_DAYS, INDIVIDUAL, SEQUENCE_INPUT, CLONE) %>%
  summarise() %>% compute_clntab() %>%
  mutate(AGE_GROUP = as.numeric(sub("-.*", "", INDIVIDUAL)),
         AGE_DAYS = age_groups[AGE_GROUP])

# Fit Zipf distributions
individuals <- tab %>% pull(INDIVIDUAL) %>% unique
n_exclude <- setNames(rep(5, length(individuals)), individuals)
n_exclude_null <- setNames(rep(0, length(individuals)), individuals)
clntab <- clntab_base %>% compute_expected_frequencies(n_exclude)
clntab_null <- clntab_base %>% compute_expected_frequencies(n_exclude_null)

# Summarise exponents
zexp <- clntab %>% group_by(AGE_DAYS, INDIVIDUAL, S) %>% summarise %>%
  mutate(S = round(S, zexp_round_ndigits))
zexp_null <- clntab_null %>% group_by(AGE_DAYS, INDIVIDUAL, S) %>% summarise %>%
  mutate(S = round(S, zexp_round_ndigits))

# Extract and save min/max exponents
zexp_min <- zexp %>% pull(S) %>% min
zexp_max <- zexp %>% pull(S) %>% max
zexp_null_min <- zexp_null %>% pull(S) %>% min
zexp_null_max <- zexp_null %>% pull(S) %>% max
savetxt(zexp_min, paste0(filename_base, "-exponent-min"))
savetxt(zexp_max, paste0(filename_base, "-exponent-max"))
savetxt(zexp_null_min, paste0(filename_base, "-exponent-null-min"))
savetxt(zexp_null_max, paste0(filename_base, "-exponent-null-max"))

# Compute actual vs expected P20
p20 <- clntab %>% filter(CLNRANK <= 20) %>% group_by(AGE_DAYS, INDIVIDUAL, S) %>% 
  summarise(Observed = sum(CLNFREQ), Expected = sum(CLNFREQ_EXP))
p20_pc <- p20 %>% mutate(Observed = Observed * 100, Expected = Expected * 100)
p20_melt <- p20_pc %>% ungroup() %>% select(-S) %>%
  melt(id.vars = "INDIVIDUAL", value.name = "P20", variable.name = "Type")
p20_obs_min <- p20_pc %>% pull(Observed) %>% min %>% round(1)
p20_obs_max <- p20_pc %>% pull(Observed) %>% max %>% round(1)
p20_exp_min <- p20_pc %>% pull(Expected) %>% min %>% round(1)
p20_exp_max <- p20_pc %>% pull(Expected) %>% max %>% round(1)
savetxt(p20_obs_min, paste0(filename_base, "-p20-obs-min"))
savetxt(p20_obs_max, paste0(filename_base, "-p20-obs-max"))
savetxt(p20_exp_min, paste0(filename_base, "-p20-exp-min"))
savetxt(p20_exp_max, paste0(filename_base, "-p20-exp-max"))

p20_cor <- cor(p20$Observed, p20$Expected) %>% round(3)

#------------------------------------------------------------------------------
# PLOT ZIFP FITS
#------------------------------------------------------------------------------

# Overlaid line plots
zplot_over <- zplot_line(clntab, "INDIVIDUAL", palette, "Individual")

# Separate point plots
zplot <- zplot_point(clntab, "INDIVIDUAL", palette, ymax = 0.1) +
  geom_text(aes(label = paste0("s = ", S), x = zexp_annot_x, y = zexp_annot_y),
            data = zexp, size = 4, family = font)
zplot_null <- zplot_point(clntab_null, "INDIVIDUAL", palette, 16, ymax = 0.1) +
  geom_text(aes(label = paste0("s = ", S), x = zexp_annot_x, y = zexp_annot_y),
            data = zexp_null, size = 4, family = font)

#------------------------------------------------------------------------------
# PLOT P20
#------------------------------------------------------------------------------

signif_stars <- function(p){
  l <- rep("n.s.", length(p))
  l <- ifelse(p <= 0.05, "*", l)
  l <- ifelse(p <= 0.01, "**", l)
  l <- ifelse(p <= 0.001, "***", l)
  return(l)
}

p20_diff <- tibble(P = wilcox.test(p20$Observed, p20$Expected)$p.value) %>%
  mutate(LABEL = signif_stars(P))

g20 <- p20_melt %>% ggplot(aes(x=Type, y=P20, fill = Type)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.7) + ylim(c(0,25)) + ylab("P20 (%)") +
  scale_fill_manual(values = c(colours[["LR"]], colours[["LG"]])) +
  geom_segment(x=1, xend = 2, y = 23, yend = 23) + 
  geom_text(x=1.5, y=23.5, size = 6, label = p20_diff$LABEL, family = font) +
  theme_classic() + theme_base + theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 0.3, unit = "cm")),
    axis.text.x = element_text(size = fontsize_base * fontscale_title,
                               margin = margin(t = 0.3, unit = "cm"))
    )

g20_cor <- p20_pc %>% ggplot() + 
  geom_point(aes(x=Observed, y=Expected, colour = INDIVIDUAL), size = 4) +
  xlab("Observed P20 (%)") + ylab("Expected P20 (%)") +
  geom_text(x=11.5, y = 6.7, label = paste0("r = ", p20_cor), size = 5,
            family = font) +
  scale_colour_manual(values = palette, name = "Individual") +
  theme_classic() + theme_base + theme(
    legend.title = element_text(margin = margin(r = 0.5, unit = "cm"))
  )

g20_plt <- gplot_grid(g20, g20_cor, nrow = 1)

#------------------------------------------------------------------------------
# SAVE PLOTS
#------------------------------------------------------------------------------

savefig(zplot_over, paste0(filename_base, "-lines"), height = 15, ratio = 1.5)
savefig(zplot, paste0(filename_base, "-fit"), height = 20, ratio = 1)
savefig(zplot_null, paste0(filename_base, "-fit-null"), height = 20, ratio = 1)
savefig(g20_plt, paste0(filename_base, "-p20-obs-exp"), height = 15, ratio = 1.8)

#------------------------------------------------------------------------------
# SAVE EXPONENT/P20 INFORMATION FOR COMPARISON WITH DIVERSITY SPECTRA
#------------------------------------------------------------------------------

# By individual
stats_indiv <- p20_pc %>% full_join(clntab_null %>% group_by(INDIVIDUAL, S) 
                                    %>% summarise, by = "INDIVIDUAL", 
                                    suffix = c("_Filtered", "_All")) %>%
  rename(P20_Observed = Observed) %>% select(-Expected)
  

# By replicate
n_exclude_rep <- setNames(rep(5,12), tab %>% pull(REPLICATE) %>% unique)
n_exclude_rep_null <- setNames(rep(0,12), tab %>% pull(REPLICATE) %>% unique)
clntab_rep_base <- tab %>% 
  group_by(INDIVIDUAL=REPLICATE, SEQUENCE_INPUT, CLONE) %>%
  summarise() %>% compute_clntab()
clntab_rep <- compute_expected_frequencies(clntab_rep_base, n_exclude_rep)
clntab_rep_null <- compute_expected_frequencies(clntab_rep_base, 
                                                n_exclude_rep_null)
stats_rep <- clntab_rep %>% group_by(INDIVIDUAL, S_Filtered = S) %>% 
  filter(CLNRANK <= 20) %>% summarise(P20_Observed = sum(CLNFREQ)) %>%
  full_join(clntab_rep_null %>% group_by(INDIVIDUAL, S) %>% summarise,
            by = "INDIVIDUAL") %>% rename(S_All = S) %>%
  ungroup() %>% rename(REPLICATE = INDIVIDUAL)

# Save
write_tsv(stats_indiv, "../_Data/changeo/pilot-clone-stats-indiv.tsv")
write_tsv(stats_rep, "../_Data/changeo/pilot-clone-stats-rep.tsv")

