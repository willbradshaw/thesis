###############################################################################
## FIGURE & TABLE                                                            ##
## Nothobranchius furzeri CH inter-sublocus dotplots                         ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")
source("aux/ggplot2.R")

# Configure input paths
ch_seqs_path_nt <- "../_Data/constant/nfu/nfu_ch_nt_exons.fasta"
ch_seqs_path_aa <- "../_Data/constant/nfu/nfu_ch_aa_exons.fasta"

# Configure output
filename_base <- "nfu-ch-aln"
filename_aa <- paste0(filename_base, "-aa")
filename_nt <- paste0(filename_base, "-nt")

#------------------------------------------------------------------------------
# IDENTIFY SUBLOCUS RANGES
#------------------------------------------------------------------------------

# Get CH sequences and exclude TM2s
ch_seqs_nt <- readDNAStringSet(ch_seqs_path_nt)
ch_seqs_aa <- readAAStringSet(ch_seqs_path_aa)

# Split by sublocus and sort ascending
ch_seqs_igh1_nt <- ch_seqs_nt[grepl("IGH1", names(ch_seqs_nt))]
ch_seqs_igh2_nt <- rev(ch_seqs_nt[grepl("IGH2", names(ch_seqs_nt))])
ch_seqs_igh1_aa <- ch_seqs_aa[grepl("IGH1", names(ch_seqs_aa))]
ch_seqs_igh2_aa <- rev(ch_seqs_aa[grepl("IGH2", names(ch_seqs_aa))])

#------------------------------------------------------------------------------
# EXTRACT SYNTENY INFORMATION AND MAKE ALIGNMENT PLOT
#------------------------------------------------------------------------------

# Pairwise identity method
pid_type <- "PID2"

# Clustering method
clust_method <- "single"

# Count D and J sequences in IGH1
n_ch_igh1_nt <- length(ch_seqs_igh1_nt)
n_ch_igh2_nt <- length(ch_seqs_igh2_nt)
n_ch_igh1_aa <- length(ch_seqs_igh1_aa)
n_ch_igh2_aa <- length(ch_seqs_igh2_aa)

# Generate pairwise alignments and compute % ID
ch_id_nt <- sapply(ch_seqs_igh1_nt, function(subject) 
  pid(pairwiseAlignment(ch_seqs_igh2_nt, subject), type = pid_type))
rownames(ch_id_nt) <- names(ch_seqs_igh2_nt)
ch_id_aa <- sapply(ch_seqs_igh1_aa, function(subject) 
  pid(pairwiseAlignment(ch_seqs_igh2_aa, subject), type = pid_type))
rownames(ch_id_aa) <- names(ch_seqs_igh2_aa)


# Convert into dataframes (separate and combined) for plotting
make_sublocus_tab <- function(id_table, seqs_igh1, seqs_igh2){
  melt(id_table, varnames = c("IGH2", "IGH1"), value.name = "ID") %>%
    mutate(IGH2 = sub("IGH2", "", IGH2), IGH1 = sub("IGH1", "", IGH1)) %>%
    mutate(IGH2 = factor(IGH2, levels = sub("IGH2", "", names(seqs_igh2))),
           IGH1 = factor(IGH1, levels = rev(sub("IGH1", "", names(seqs_igh1)))),
           ISO1 = sub("-.*", "", IGH1), 
           ISO2 = sub("-.*", "", IGH2),
           ISO = ifelse(ISO1 != ISO2, "GR", paste0("C", ISO1)),
           ISO = factor(ISO, levels = names(colours))
    ) %>% as.tibble %>% select(-ISO1, -ISO2)
}
ch_tab_nt <- make_sublocus_tab(ch_id_nt, ch_seqs_igh1_nt, ch_seqs_igh2_nt)
ch_tab_aa <- make_sublocus_tab(ch_id_aa, ch_seqs_igh1_aa, ch_seqs_igh2_aa)
ch_tab_comb <- bind_rows(ch_tab_nt %>% mutate(TYPE = "Nucleotide sequence identity"),
                         ch_tab_aa %>% mutate(TYPE = "Amino-acid sequence identity"))

# Define colour scheme and theme (same for all plots)
label_colours <- rownames(ch_id_nt) %>% .[!grepl("TM2", .)] %>% 
  sub("IGH\\d(.*)-.*", "C\\1", .) %>% colours[.] %>% as.character
theme_hm <- theme_classic() + theme_base + theme(
  axis.text.x = element_text(angle = 90, colour = label_colours,
                             vjust = 0.5, hjust = 1),
  axis.text.y = element_text(colour = rev(label_colours),
                             vjust = 0.5, hjust = 1)
)

# Cut tables at minimum visualised identity and exclude TM2 exons
range_min <- 60
ch_tab_cut_nt <- ch_tab_nt %>% mutate(ID = pmax(range_min, ID)) %>%
  filter(!grepl("TM2", IGH1), !grepl("TM2", IGH2))
ch_tab_cut_aa <- ch_tab_aa %>% mutate(ID = pmax(range_min, ID)) %>%
  filter(!grepl("TM2", IGH1), !grepl("TM2", IGH2))
ch_tab_cut_comb <- ch_tab_comb %>% mutate(ID = pmax(range_min, ID)) %>%
  filter(!grepl("TM2", IGH1), !grepl("TM2", IGH2))

ch_heatmap <- function(tab){
  ggplot(tab) + geom_tile(aes(x=IGH2, y=IGH1, fill=ISO, alpha=ID)) +
  scale_fill_manual(values = unlist(colours), name="Isotype", labels = c("IGHD", "IGHM", "Unmatched")) +
  scale_alpha_continuous(name="% ID", limits=c(range_min,100)) + 
  coord_fixed() + theme_hm
}

g_ch_nt <- ch_heatmap(ch_tab_cut_nt)
g_ch_aa <- ch_heatmap(ch_tab_cut_aa)
g_ch_comb <- ch_heatmap(ch_tab_cut_comb) + facet_grid(.~TYPE)

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plot_height <- 20
plot_ratio <- 1/1.1

savefig(g_ch_nt, filename_nt, height = plot_height, ratio = plot_ratio)
savefig(g_ch_aa, filename_aa, height = plot_height, ratio = plot_ratio)
savefig(g_ch_comb, filename_base, height = plot_height, ratio = plot_ratio * 1.65)

#------------------------------------------------------------------------------
# GENERATE LATEX TABLE
#------------------------------------------------------------------------------

defactor <- function(tab){
  tab %>% mutate(IGH1 = as.character(IGH1), 
                 IGH2 = as.character(IGH2),
                 ISO = as.character(ISO))
}
ch_retab_nt <- defactor(ch_tab_nt) %>% rename(NT = ID) %>% select(-ISO)
ch_retab_aa <- defactor(ch_tab_aa) %>% rename(AA = ID) %>% select(-ISO)

ch_out_tab <- full_join(ch_retab_nt, ch_retab_aa, by=c("IGH1", "IGH2")) %>%
  filter(IGH1 == IGH2) %>%
  mutate(Isotype = sub("(.*)-.*", "\\1", IGH1),
         Exon = sub(".*-(.*)", "\\1", IGH1),
         NT = round(NT, 2),
         AA = round(AA, 2)) %>%
  select(Isotype, Exon, NT, AA)

savetab(ch_out_tab, filename_base)

#------------------------------------------------------------------------------
# CALCULATE AVERAGE IDENTITIES
#------------------------------------------------------------------------------
# TODO: Sweave this?

ch_avg_tab <- ch_out_tab %>% group_by(Isotype) %>%
  summarise(NT = mean(NT), AA = mean(AA))