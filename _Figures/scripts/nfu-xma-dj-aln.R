###############################################################################
## FIGURE                                                                    ##
## N. furzeri vs X. maculatus D/J inter-sublocus dotplots                    ##
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
dh_seqs_path_nfu <- "../_Data/segments/nfu/nfu_dh_nt.fasta"
jh_seqs_path_nfu <- "../_Data/segments/nfu/nfu_jh_nt.fasta"
dh_seqs_path_xma <- "../_Data/segments/xma/xma_dh_nt.fasta"
jh_seqs_path_xma <- "../_Data/segments/xma/xma_jh_nt.fasta"

# Configure output
filename <- "nfu-xma-dj-aln"

#------------------------------------------------------------------------------
# IDENTIFY SUBLOCUS RANGES
#------------------------------------------------------------------------------

# Get D/J sequences
get_sorted_segments <- function(path){
  readDNAStringSet(path) #%>% .[match(sort(names(.)), names(.))]
}
dh_seqs_nfu <- get_sorted_segments(dh_seqs_path_nfu)
jh_seqs_nfu <- get_sorted_segments(jh_seqs_path_nfu)
dh_seqs_xma <- get_sorted_segments(dh_seqs_path_xma)
jh_seqs_xma <- get_sorted_segments(jh_seqs_path_xma)

#------------------------------------------------------------------------------
# EXTRACT SYNTENY INFORMATION AND MAKE ALIGNMENT PLOT
#------------------------------------------------------------------------------

# Pairwise identity method
pid_type <- "PID2"

# Clustering method
clust_method <- "single"

# Count D and J sequences
n_dh_nfu <- length(dh_seqs_nfu)
n_dh_xma <- length(dh_seqs_xma)
n_jh_nfu <- length(jh_seqs_nfu)
n_jh_xma <- length(jh_seqs_xma)


# Generate pairwise alignments and compute % ID
dh_id <- sapply(dh_seqs_nfu, function(subject) 
  pid(pairwiseAlignment(dh_seqs_xma, subject), type = pid_type))
rownames(dh_id) <- names(dh_seqs_xma)
jh_id <- sapply(jh_seqs_nfu, function(subject) 
  pid(pairwiseAlignment(jh_seqs_xma, subject), type = pid_type))
rownames(jh_id) <- names(jh_seqs_xma)

# Convert into dataframe for plotting
dh_tab <- melt(dh_id, varnames = c("XMA", "NFU"), value.name = "ID") %>%
  mutate(XMA = factor(XMA, levels = names(dh_seqs_xma)),
         NFU = factor(NFU, levels = rev(names(dh_seqs_nfu))))
jh_tab <- melt(jh_id, varnames = c("XMA", "NFU"), value.name = "ID") %>%
  mutate(XMA = factor(XMA, levels = names(jh_seqs_xma)),
         NFU = factor(NFU, levels = rev(names(jh_seqs_nfu))))

# Determine positions of sublocus and Z/M cutoffs
n_dh_igh2_nfu <- sum(grepl("IGH2", names(dh_seqs_nfu)))
n_jh_igh2_nfu <- sum(grepl("IGH2", names(jh_seqs_nfu)))
n_dh_z_xma <- sum(grepl("Z", names(dh_seqs_xma)))
n_jh_z_xma <- sum(grepl("Z", names(jh_seqs_xma)))


# Make identity heatmaps
theme_hm <- theme_classic() + theme_base +
  theme(axis.text.x = element_text(angle=90))

range_min <- 60
dh_tab_cut <- dh_tab %>% mutate(ID = pmax(range_min, ID))
jh_tab_cut <- jh_tab %>% mutate(ID = pmax(range_min, ID))

g_dh <- ggplot(dh_tab_cut) + 
  geom_tile(aes(x=XMA, y=NFU, fill=ID)) + 
  scale_fill_gradient(low="white", high=colours["DH"],
                      name="% ID", limits=c(range_min,100)) + 
  geom_hline(yintercept = n_dh_igh2_nfu + 0.5) +
  geom_vline(xintercept = n_dh_z_xma + 0.5) +
  coord_fixed() + theme_hm
g_jh <- ggplot(jh_tab_cut) + 
  geom_tile(aes(x=XMA, y=NFU, fill=ID)) + 
  scale_fill_gradient(low="white", high=colours["JH"],
                      name="% ID", limits=c(range_min,100)) + 
  geom_hline(yintercept = n_jh_igh2_nfu + 0.5) +
  geom_vline(xintercept = n_jh_z_xma + 0.5) +
  coord_fixed() + theme_hm

width_adjust <- 18
plt <- plot_grid(g_dh, g_jh, ncol=2,
                 rel_widths = c(n_dh_xma+width_adjust, n_jh_xma+width_adjust),
                 labels="AUTO", label_fontfamily = titlefont,
                 label_fontface = "plain",
                 label_size = fontsize_base * fontscale_label
)
plt

#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

plot_height <- 20
plot_ratio <- 1.3

savefig(plot = plt, filename = filename, device = "svg", 
        height = plot_height, ratio = plot_ratio)
savefig(plot = plt, filename = filename, device = "png", 
        height = plot_height, ratio = plot_ratio)


