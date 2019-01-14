###############################################################################
## FIGURE & TABLE                                                            ##
## Inter-IGHZ sequence similarity between P. playfairii and other species    ##
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

# Configure input paths for CZ-bearing species
ch_paths_nt <- character(0)
ch_paths_aa <- character(0)
species <- c("xma", "pre", "pfo", "fhe", "cva", "kma", "ppl", "cto")
ch_paths_nt <- sapply(species, function(s) 
  paste0("../_Data/constant/", s, "/", s, "_ch_nt_4trim.fasta"))
ch_paths_aa <- sapply(species, function(s) 
  paste0("../_Data/constant/", s, "/", s, "_ch_aa.fasta"))
ch_paths_nt <- gsub("xma_ch", "xma_new_ch", ch_paths_nt)
ch_paths_aa <- gsub("xma_ch", "xma_new_ch", ch_paths_aa)

# Configure output
filename_base <- "ppl-cz-aln"
filename_aa <- paste0(filename_base, "-aa")
filename_nt <- paste0(filename_base, "-nt")

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

# Determine species identity of sequences
ch_nt_list <- lapply(ch_paths_nt, readDNAStringSet)
ch_aa_list <- lapply(ch_paths_aa, readAAStringSet)
species_nt <- unlist(mapply(rep, names(ch_nt_list), sapply(ch_nt_list, length)))
species_aa <- unlist(mapply(rep, names(ch_aa_list), sapply(ch_aa_list, length)))

# Import and annotate sequences
ch_nt <- readDNAStringSet(ch_paths_nt)
names(ch_nt) <- paste0(species_nt, "_", names(ch_nt))
ch_aa <- readAAStringSet(ch_paths_aa)
names(ch_aa) <- paste0(species_aa, "_", names(ch_aa))

# Subset to CZ
cz_pattern <- "IGH.?Z.?-[1-4]$"
cz_nt <- ch_nt[grepl(cz_pattern, names(ch_nt))]
cz_aa <- ch_aa[grepl(cz_pattern, names(ch_aa))]

# Separate ppl sequences
cz_nt_ppl <- cz_nt[grepl("ppl", names(cz_nt))]
cz_aa_ppl <- cz_aa[grepl("ppl", names(cz_aa))]
cz_nt_other <- cz_nt[!grepl("ppl", names(cz_nt))]
cz_aa_other <- cz_aa[!grepl("ppl", names(cz_aa))]

#------------------------------------------------------------------------------
# ALIGN EXON PAIRS
#------------------------------------------------------------------------------

# Compute pairwise alignments and return scores
cz_id_nt <- sapply(cz_nt_ppl, function(subject) 
  pairwiseAlignment(cz_nt_other, subject, scoreOnly = TRUE))
rownames(cz_id_nt) <- names(cz_nt_other)
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
    as.tibble %>% select(-PPL, -OTHER)
}
cz_tab_nt <- make_tab(cz_id_nt)
cz_tab_aa <- make_tab(cz_id_aa)
cz_tab_comb <- bind_rows(cz_tab_nt %>% mutate(TYPE = "NT"),
                         cz_tab_aa %>% mutate(TYPE = "AA"))

#------------------------------------------------------------------------------
# PLOT SCORE DISTRIBUTIONS
#------------------------------------------------------------------------------

cz_score_box <- function(tab){
  ggplot(tab) + 
    geom_boxplot(aes(x = factor(EXON_PPL), y = SCORE, colour = LINEAGE_OTHER),
                 width = 0.4, position = position_dodge(0.5)) +
    xlab("CZ exon") + ylab("Alignment score") + 
    scale_colour_discrete(name = "IGHZ lineage") +
    theme_base + theme(legend.justification = "centre")
}
g_cz_nt <- cz_score_box(cz_tab_nt)
g_cz_aa <- cz_score_box(cz_tab_aa)
g_cz_comb <- cz_score_box(cz_tab_comb) + facet_wrap(~TYPE, ncol=1)

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plot_height <- 20
plot_ratio <- 1/1.3

savefig(g_cz_nt, filename_nt, height = plot_height, ratio = plot_ratio)
savefig(g_cz_aa, filename_aa, height = plot_height, ratio = plot_ratio)
savefig(g_cz_comb, filename_base, height = plot_height, ratio = plot_ratio * 1.65)