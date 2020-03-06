###############################################################################
## FIGURE & TABLE                                                            ##
## Inter-IGHZ sequence similarity in X. maculatus IGH locus                  ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# IDENTIFY SUBLOCUS RANGES
#------------------------------------------------------------------------------

# Get CH sequences and exclude TM2s
ch_seqs_nt <- readDNAStringSet(inpath_nt)
ch_seqs_aa <- readAAStringSet(inpath_aa)

# Extract IGHZ1 and IGHZ2 exons
cz1_seqs_nt <- ch_seqs_nt[grepl("IGHZ1", names(ch_seqs_nt))]
cz1_seqs_aa <- ch_seqs_aa[grepl("IGHZ1", names(ch_seqs_aa))]
cz2_seqs_nt <- ch_seqs_nt[grepl("IGHZ2", names(ch_seqs_nt))]
cz2_seqs_aa <- ch_seqs_aa[grepl("IGHZ2", names(ch_seqs_aa))]

#------------------------------------------------------------------------------
# ALIGN EXON PAIRS
#------------------------------------------------------------------------------

# Count D and J sequences in IGH1
n_cz1_nt <- length(cz1_seqs_nt)
n_cz1_aa <- length(cz1_seqs_aa)
n_cz2_nt <- length(cz2_seqs_nt)
n_cz2_aa <- length(cz2_seqs_aa)

# Generate pairwise alignments and compute % ID
cz_id_nt <- sapply(cz1_seqs_nt, function(subject) 
  pid(pairwiseAlignment(cz2_seqs_nt, subject), type = pid_type))
rownames(cz_id_nt) <- names(cz2_seqs_nt)
cz_id_aa <- sapply(cz1_seqs_aa, function(subject) 
  pid(pairwiseAlignment(cz2_seqs_aa, subject), type = pid_type))
rownames(cz_id_aa) <- names(cz2_seqs_aa)

# Convert into dataframes (separate and combined) for plotting
make_sublocus_tab <- function(id_table, seqs_ighz1, seqs_ighz2){
  melt(id_table, varnames = c("IGHZ2", "IGHZ1"), value.name = "ID") %>%
    mutate(IGHZ2 = sub("IGHZ2-", "", IGHZ2), IGHZ1 = sub("IGHZ1-", "", IGHZ1)) %>%
    mutate(IGHZ2 = factor(IGHZ2, levels = sub("IGHZ2-", "", names(seqs_ighz2))),
           IGHZ1 = factor(IGHZ1, levels = rev(sub("IGHZ1-", "", names(seqs_ighz1)))),
           ISO1 = "Z", 
           ISO2 = "Z",
           ISO = ifelse(ISO1 != ISO2, "GR", paste0("C", ISO1)),
           ISO = factor(ISO, levels = names(colours))
    ) %>% as.tibble %>% select(-ISO1, -ISO2)
}
cz_tab_nt <- make_sublocus_tab(cz_id_nt, cz1_seqs_nt, cz2_seqs_nt)
cz_tab_aa <- make_sublocus_tab(cz_id_aa, cz1_seqs_aa, cz2_seqs_aa)
cz_tab_comb <- bind_rows(cz_tab_nt %>% mutate(TYPE = "Nucleotide sequence identity"),
                         cz_tab_aa %>% mutate(TYPE = "Amino-acid sequence identity"))


# Define colour scheme and theme (same for all plots)
label_colours <- rownames(cz_id_nt) %>% .[!grepl("TM[12]", .)] %>% 
  .[!grepl("S$", .)] %>% 
  sub("IGH\\d?(.)\\d?-.*", "C\\1", .) %>% colours[.] %>% as.character
theme_hm <- theme_classic() + theme_base + theme(
  axis.text.x = element_text(colour = label_colours,
                             vjust = 0.5, hjust = 0.5),
  axis.text.y = element_text(colour = rev(label_colours),
                             vjust = 0.5, hjust = 1)
)

# Cut tables at minimum visualised identity and exclude TM2 exons
range_min <- 20
range_max <- 70
cz_tab_cut_nt <- cz_tab_nt %>% mutate(ID = pmax(range_min, ID)) %>%
  filter(!grepl("TM[12]", IGHZ1), !grepl("TM[12]", IGHZ2),
         !grepl("S$", IGHZ1), !grepl("S$", IGHZ2))
cz_tab_cut_aa <- cz_tab_aa %>% mutate(ID = pmax(range_min, ID)) %>%
  filter(!grepl("TM[12]", IGHZ1), !grepl("TM[12]", IGHZ2),
         !grepl("S$", IGHZ1), !grepl("S$", IGHZ2))
cz_tab_cut_comb <- cz_tab_comb %>% mutate(ID = pmax(range_min, ID)) %>%
  filter(!grepl("TM[12]", IGHZ1), !grepl("TM[12]", IGHZ2),
         !grepl("S$", IGHZ1), !grepl("S$", IGHZ2))

ch_heatmap <- function(tab){
  ggplot(tab) + geom_tile(aes(x=IGHZ2, y=IGHZ1, fill=ISO, alpha=ID)) +
  scale_fill_manual(values = unlist(colours), name="Isotype", labels = c("IGHZ")) +
  scale_alpha_continuous(name="% ID", limits=c(range_min,range_max)) + 
  coord_fixed() + theme_hm
}

g_cz_nt <- ch_heatmap(cz_tab_cut_nt)
g_cz_aa <- ch_heatmap(cz_tab_cut_aa)
g_cz_comb <- ch_heatmap(cz_tab_cut_comb) + facet_grid(.~TYPE)

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plot_height <- 20
plot_ratio <- 1/1.1

ggsave(plot = g_cz_comb, filename = outpath_map, device = "svg", units = "cm",
       height=plot_height, width=plot_height * (plot_ratio*1.65))

#------------------------------------------------------------------------------
# GENERATE LATEX TABLE
#------------------------------------------------------------------------------

defactor <- function(tab){
  tab %>% mutate(IGHZ1 = as.character(IGHZ1),
                 IGHZ2 = as.character(IGHZ2),
                 ISO = as.character(ISO))
}
cz_retab_nt <- defactor(cz_tab_nt) %>% rename(NT = ID) %>% select(-ISO)
cz_retab_aa <- defactor(cz_tab_aa) %>% rename(AA = ID) %>% select(-ISO)

cz_out_tab <- full_join(cz_retab_nt, cz_retab_aa, by=c("IGHZ1", "IGHZ2")) %>%
  filter(IGHZ1 == IGHZ2, 
         !IGHZ1 %in% c("S", "TM1", "TM2")) %>%
  mutate(Isotype = "Z",
         Exon = IGHZ1,
         NT = round(NT, 2),
         AA = round(AA, 2)) %>%
  select(Isotype, Exon, NT, AA)

savetab(cz_out_tab, outpath_tab)

#------------------------------------------------------------------------------
# CALCULATE AVERAGE IDENTITIES
#------------------------------------------------------------------------------

cz_avg_tab <- cz_out_tab %>% group_by(Isotype) %>%
  filter(! Exon %in% c("S", "TM1", "TM2")) %>%
  summarise(NT = mean(NT), AA = mean(AA))
