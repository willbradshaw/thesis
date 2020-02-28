###############################################################################
## FIGURE & TABLE                                                            ##
## Inter-IGHZ sequence similarity between P. playfairii and other species    ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Specify lineages
ighz_lineages <- list(
  "pre_IGHZ1" = "A",
  "pfo_IGHZ1" = "A",
  "xma_IGHZ1" = "A",
  "fhe_IGHZ1" = "A",
  "kma_IGHZ1" = "A",
  "cto_IGHZ3" = "B",
  "fhe_IGHZ2" = "C",
  "cva_IGHZ"  = "C",
  "pfo_IGHZ3" = "C",
  "xma_IGHZ2" = "C",
  "kma_IGHZ2" = "C",
  "pfo_IGHZ2" = "B",
  "pre_IGHZ2" = "B",
  "cto_IGHZ1" = "A",
  "cto_IGHZ2" = "B",
  "cto_IGHZ4" = "B"
)

#------------------------------------------------------------------------------
# IDENTIFY SUBLOCUS RANGES
#------------------------------------------------------------------------------
names(inpaths) <- species

# Determine species identity of sequences
ch_aa_list <- lapply(inpaths, readAAStringSet)
species_aa <- unlist(mapply(rep, names(ch_aa_list), sapply(ch_aa_list, length)))

# Import and annotate sequences
ch_aa <- readAAStringSet(unlist(inpaths))
names(ch_aa) <- paste0(species_aa, "_", names(ch_aa))

# Subset to CZ
cz_pattern <- "IGH.?Z.?-[1-4]$"
cz_aa <- ch_aa[grepl(cz_pattern, names(ch_aa))]

# Separate ppl sequences
cz_aa_ppl <- cz_aa[grepl("ppl", names(cz_aa))]
cz_aa_other <- cz_aa[!grepl("ppl", names(cz_aa))]

#------------------------------------------------------------------------------
# ALIGN EXON PAIRS
#------------------------------------------------------------------------------
# Compute pairwise alignments and return scores
cz_id_aa <- sapply(cz_aa_ppl, function(subject) 
  pairwiseAlignment(cz_aa_other, subject, scoreOnly = TRUE))
rownames(cz_id_aa) <- names(cz_aa_other)

# Convert into dataframes (separate and combined) for plotting
make_tab <- function(score_table){
  melt(score_table, varnames = c("OTHER", "PPL"), value.name = "SCORE") %>%
    mutate(EXON_PPL = as.integer(sub("ppl_IGHZ-", "", PPL)),
           SPECIES_OTHER = sub("_.*", "", OTHER),
           EXON_OTHER = as.integer(sub(".*IGHZ\\d?-(\\d*).*", "\\1", OTHER)),
           LOCUS_OTHER = sub(".*IGHZ(.*)-\\d.*", "\\1", OTHER),
           LINEAGE_OTHER = unlist(ighz_lineages[paste0(SPECIES_OTHER, "_", 
                                                 "IGHZ", LOCUS_OTHER)])
    ) %>%
    filter(EXON_PPL == EXON_OTHER) %>%
    as_tibble %>% select(-PPL, -OTHER)
}
cz_tab_aa <- make_tab(cz_id_aa)

#------------------------------------------------------------------------------
# TEST FOR SIGNIFICANCE
#------------------------------------------------------------------------------
# Prepare significance annotations
signif_stars <- function(p){
  l <- rep("n.s.", length(p))
  l <- ifelse(p <= 0.05, "*", l)
  l <- ifelse(p <= 0.01, "**", l)
  l <- ifelse(p <= 0.001, "***", l)
  return(l)
}
annotate_sigtable <- function(ztab, scale = 8, width = 0.16){
  ptab <- lapply(seq(3), function(x) lapply(seq(4), function(n) 
    ztab %>% filter(EXON_PPL == n, LINEAGE_OTHER %in% LETTERS[1:3][-x]) %>%
      wilcox.test(formula = SCORE ~ LINEAGE_OTHER, data = .) %>%
      (function(x) x$p.value) %>%
      tibble(EXON_PPL = n, P = ., LIN1 = LETTERS[1:3][-x][1],
             LIN2 = LETTERS[1:3][-x][2])) %>% bind_rows) %>% bind_rows %>%
    filter(P <= significance_level)
  ytab <- ztab %>% group_by(EXON_PPL, LINEAGE_OTHER) %>%
    summarise(SCORE = max(SCORE))
  qtab <- ytab %>% rename(LIN1 = LINEAGE_OTHER, Y1 = SCORE) %>%
    inner_join(ptab, by = c("EXON_PPL", "LIN1")) %>%
    inner_join(ytab %>% rename(LIN2 = LINEAGE_OTHER, Y2 = SCORE),
             by = c("EXON_PPL", "LIN2")) %>%
    group_by(EXON_PPL) %>% arrange(LIN1, LIN2) %>%
  mutate(Y = pmax(Y1, Y2) + scale * (1 + row_number() / 2),
         Y_LAB = Y + (scale * 0.3),
         X1 = EXON_PPL + ifelse(LIN1 == "A", -width, 0),
         X2 = EXON_PPL + ifelse(LIN2 == "C", width, 0),
         X_AVG = (X1 + X2)/2,
         LABEL = signif_stars(P))
  # select(-Y1, -Y2)
  return(qtab)
}

p_aa <- annotate_sigtable(cz_tab_aa)

         
#------------------------------------------------------------------------------
# PLOT SCORE DISTRIBUTIONS
#------------------------------------------------------------------------------

cz_score_box <- function(ztab, ptab){
  ggplot(ztab) + 
    geom_boxplot(aes(x = factor(EXON_PPL), y = SCORE, colour = LINEAGE_OTHER),
                 width = 0.4, position = position_dodge(0.5)) +
    geom_segment(aes(x=X1, xend=X2, y=Y, yend=Y), data = ptab) +
    geom_text(aes(x=X_AVG, y=Y_LAB, label=LABEL), data = ptab, size = 6) +
    xlab("CZ exon") + ylab("Alignment score") + 
    scale_colour_discrete(name = "IGHZ lineage") +
    theme_classic() + theme_base + theme(legend.justification = "centre")
}
g_cz_aa <- cz_score_box(cz_tab_aa, p_aa)
#g_cz_comb <- cz_score_box(cz_tab_comb) + facet_wrap(~TYPE, ncol=1)

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plot_height <- 12
plot_width <- 18

ggsave(plot = g_cz_aa, filename = outpath, device = "svg", units = "cm",
       height = plot_height, width = plot_width)
