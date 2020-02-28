###############################################################################
## FIGURE                                                                    ##
## Nothobranchius furzeri D/J inter-sublocus dotplots                        ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# IDENTIFY SUBLOCUS RANGES
#------------------------------------------------------------------------------

# Get D/J sequences
dh_seqs <- readDNAStringSet(inpath_dh)
jh_seqs <- readDNAStringSet(inpath_jh)

# Split by sublocus and sort ascending
dh_seqs_igh1 <- dh_seqs[grepl("IGH1", names(dh_seqs))]
dh_seqs_igh1 <- dh_seqs_igh1[match(names(dh_seqs_igh1), sort(names(dh_seqs_igh1)))]
dh_seqs_igh2 <- dh_seqs[grepl("IGH2", names(dh_seqs))]
dh_seqs_igh2 <- dh_seqs_igh2[match(names(dh_seqs_igh2), sort(names(dh_seqs_igh2)))]
jh_seqs_igh1 <- jh_seqs[grepl("IGH1", names(jh_seqs))]
jh_seqs_igh1 <- jh_seqs_igh1[match(names(jh_seqs_igh1), sort(names(jh_seqs_igh1)))]
jh_seqs_igh2 <- jh_seqs[grepl("IGH2", names(jh_seqs))]
jh_seqs_igh2 <- jh_seqs_igh2[match(names(jh_seqs_igh2), sort(names(jh_seqs_igh2)))]


#------------------------------------------------------------------------------
# EXTRACT SYNTENY INFORMATION AND MAKE ALIGNMENT PLOT
#------------------------------------------------------------------------------

# Pairwise identity method
pid_type <- "PID2"

# Clustering method
clust_method <- "single"

# Count D and J sequences in IGH1
n_dh_igh1 <- length(dh_seqs_igh1)
n_dh_igh2 <- length(dh_seqs_igh2)
n_jh_igh1 <- length(jh_seqs_igh1)
n_jh_igh2 <- length(jh_seqs_igh2)


# Generate pairwise alignments and compute % ID
dh_id <- sapply(dh_seqs_igh1, function(subject) 
  pid(pairwiseAlignment(dh_seqs_igh2, subject), type = pid_type))
rownames(dh_id) <- names(dh_seqs_igh2)
jh_id <- sapply(jh_seqs_igh1, function(subject) 
  pid(pairwiseAlignment(jh_seqs_igh2, subject), type = pid_type))
rownames(jh_id) <- names(jh_seqs_igh2)

# Convert into dataframe for plotting
dh_tab <- melt(dh_id, varnames = c("IGH2D", "IGH1D"), value.name = "ID") %>%
  mutate(IGH2D = sub("IGH2D", "", IGH2D), IGH1D = sub("IGH1D", "", IGH1D))
jh_tab <- melt(jh_id, varnames = c("IGH2J", "IGH1J"), value.name = "ID") %>%
  mutate(IGH2J = sub("IGH2J", "", IGH2J), IGH1J = sub("IGH1J", "", IGH1J))

# Make identity heatmaps
theme_hm <- theme_classic() + theme_base

range_min <- 60
dh_tab_cut <- dh_tab %>% mutate(ID = pmax(range_min, ID))
jh_tab_cut <- jh_tab %>% mutate(ID = pmax(range_min, ID))

g_dh <- ggplot(dh_tab_cut) + 
  geom_tile(aes(x=IGH2D, y=IGH1D, fill=ID)) + 
  scale_fill_gradient(low="white", high=colours["DH"],
                      name="% ID", limits=c(range_min,100)) + 
  coord_fixed() + theme_hm
g_jh <- ggplot(jh_tab_cut) + 
  geom_tile(aes(x=IGH2J, y=IGH1J, fill=ID)) + 
  scale_fill_gradient(low="white", high=colours["JH"],
                      name="% ID", limits=c(range_min,100)) + 
  coord_fixed() + theme_hm

plt <- plot_grid(g_dh, g_jh, ncol=2,
                 rel_widths = c(n_dh_igh2+0.4, n_jh_igh2+0.4),
                 labels="AUTO", label_fontfamily = titlefont,
                 label_fontface = "plain",
                 label_size = fontsize_base * fontscale_label
)

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plot_height <- 15
plot_ratio <- 1.3

ggsave(plot =plt, filename = outpath, device = "svg", height=plot_height,
       width = plot_height*plot_ratio, units = "cm")
