###############################################################################
## FIGURE                                                                    ##
## Xiphophorus maculatus RSS composition plots                               ##
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

# Configure input paths
vh_rss_path_xma <- "../_Data/rss/xma/xma_new_vh_rss.fasta"
dh_rss_path_xma <- "../_Data/rss/xma/xma_new_dh_rss.fasta"
jh_rss_path_xma <- "../_Data/rss/xma/xma_new_jh_rss.fasta"

# Configure output
filename_xma_all <- "xma-new-rss-seqlogo-all"
filename_xma_sep <- "xma-new-rss-seqlogo-sep"

#------------------------------------------------------------------------------
# IMPORT RSS SEQUENCES AND EXTRACT HEPTAMER/NONAMER/SPACER SEQUENCES
#------------------------------------------------------------------------------

# Get D/J sequences
vh_rss_xma <- readDNAStringSet(vh_rss_path_xma) %>%
  .[width(.)>0]
dh_rss_xma <- readDNAStringSet(dh_rss_path_xma)
jh_rss_xma <- readDNAStringSet(jh_rss_path_xma)

# Group by length
rss_xma <- list(
  vh = vh_rss_xma,
  dh = dh_rss_xma,
  jh = jh_rss_xma,
  long = c(vh_rss_xma, jh_rss_xma),
  short = dh_rss_xma,
  all = c(vh_rss_xma, dh_rss_xma, jh_rss_xma)
)

# Extract subsequences

heptamers_xma <- lapply(rss_xma, function(seqs) subseq(seqs, 1, 7))
nonamers_xma <- lapply(rss_xma, function(seqs) subseq(seqs, -9, -1))
spacers_xma <- lapply(rss_xma, function(seqs) subseq(seqs, 8, -10))

#------------------------------------------------------------------------------
# COMPUTE PWMS AND LENGTH DISTRIBUTIONS
#------------------------------------------------------------------------------

# Character vectors
charHeptamer_xma <- lapply(heptamers_xma, function(seqs) as.character(seqs))
charNonamer_xma <- lapply(nonamers_xma, function(seqs) as.character(seqs))

# Spacer length tables
getLengthTable <- function(stringSet, type, length){
  tibble(name = names(stringSet), width = width(stringSet),
                seq = as.character(stringSet), type=type, length=length)
}
lengths <- list(vh="long", dh="short", jh="long", long="long", short="short", all="all")
spacerTable_xma <- bind_rows(
  getLengthTable(spacers_xma$vh, "VH", "long"),
  getLengthTable(spacers_xma$dh, "DH", "short"),
  getLengthTable(spacers_xma$jh, "JH", "long")
) %>% mutate(type = factor(type, levels = names(colours)))

#------------------------------------------------------------------------------
# MAKE GGPLOTS
#------------------------------------------------------------------------------

# Define theme
theme_logo <- theme_classic() + theme_base + theme(legend.position = "none")

# Sequence logos
logo_method <- "bits" # Or "probability"
logo_seq_type <- "dna"
logo_font <- "roboto_regular" # see list_fonts() for options
logo_col_scheme <- "nucleotide" # or nucleotide2, or a custom scheme

gglogo <- function(data){
  geom_logo(data = data, method = logo_method, seq_type = logo_seq_type, 
            font = logo_font, col_scheme = logo_col_scheme)
}
ggHeptamer_xma <- lapply(charHeptamer_xma, function(seqs) ggplot() + gglogo(seqs) +
                       xlab("Position") + theme_logo)
ggNonamer_xma <- lapply(charNonamer_xma, function(seqs) ggplot() + gglogo(seqs) +
                      xlab("Position") + theme_logo)


# Length distributions
tick_colours = rep(c("black", "red", "black"),2)

ggSpacer_xma <- ggplot(spacerTable_xma) + 
  geom_freqpoly(aes(x=width, colour=type), binwidth=1, size=2) +
  xlab("Length (nt)") + ylab("Frequency") +
  scale_colour_manual(values = as.character(colours), name = "Segment type") +
  scale_x_continuous(breaks = c(10, 12, 15, 20, 23, 25), limits = c(9,27)) +
  geom_vline(xintercept=12, colour = "red", linetype = "dashed") + 
  geom_vline(xintercept=23, colour = "red", linetype = "dashed") +
  theme_logo + 
  theme(axis.text.x = element_text(colour = tick_colours),
        axis.ticks.x = element_line(colour = tick_colours))

#------------------------------------------------------------------------------
# THREE-PART ALL-SEQUENCE PLOT
#------------------------------------------------------------------------------
annot_coords_xma <- list(
  VH = list(x = 19.8, y = 14),
  DH = list(x = 14.2, y = 23),
  JH = list(x = 26, y = 12)
)

annotate_spacer_plot <- function(spacer_plot, annot_coords){
  spacer_plot + 
  annotate("text", label="VH", colour = colours[["VH"]],
           x = annot_coords[["VH"]][["x"]], 
           y = annot_coords[["VH"]][["y"]], 
           size = fontsize_base/fontscale_title) +
  annotate("text", label="DH", colour = colours[["DH"]],
           x = annot_coords[["DH"]][["x"]], 
           y = annot_coords[["DH"]][["y"]], 
           size = fontsize_base/fontscale_title) +
  annotate("text", label="JH", colour = colours[["JH"]],
           x = annot_coords[["JH"]][["x"]], 
           y = annot_coords[["JH"]][["y"]], 
           size = fontsize_base/fontscale_title)
}

spacer_all_annotated_xma <- annotate_spacer_plot(ggSpacer_xma, annot_coords_xma)

plt_xma_all <- plot_grid(ggHeptamer_xma$all + ggtitle("Heptamer composition"),
             spacer_all_annotated_xma + ggtitle("Spacer length"),
             ggNonamer_xma$all + ggtitle("Nonamer composition"),
             ncol = 3, labels="AUTO", label_fontfamily = titlefont,
             label_fontface = "plain",
             label_size = fontsize_base * fontscale_label)

#------------------------------------------------------------------------------
# SPLIT-BY-SEGMENT-TYPE SEQLOGO PLOT
#------------------------------------------------------------------------------

plt_xma_sep <- plot_grid(ggHeptamer_xma$vh + ggtitle("VH heptamers"),
                         ggNonamer_xma$vh + ggtitle("VH nonamers"),
                         ggHeptamer_xma$dh + ggtitle("DH heptamers"),
                         ggNonamer_xma$dh + ggtitle("DH nonamers"),
                         ggHeptamer_xma$jh + ggtitle("JH heptamers"),
                         ggNonamer_xma$jh + ggtitle("JH nonamers"),
                         ncol = 2, labels="AUTO", label_fontfamily = titlefont,
                         label_fontface = "plain",
                         label_size = fontsize_base * fontscale_label)

#------------------------------------------------------------------------------
# SAVE PLOTS
#------------------------------------------------------------------------------

plot_all_height <- 10
plot_all_ratio <- 3 # Width vs height

savefig(plot = plt_xma_all, filename = filename_xma_all, device = "svg", 
        height = plot_all_height, ratio = plot_all_ratio)
savefig(plot = plt_xma_all, filename = filename_xma_all, device = "png", 
        height = plot_all_height, ratio = plot_all_ratio)

plot_sep_height <- 30
plot_sep_ratio <- 2/3 # Width vs height

savefig(plot = plt_xma_sep, filename = filename_xma_sep, device = "svg", 
        height = plot_sep_height, ratio = plot_sep_ratio)
savefig(plot = plt_xma_sep, filename = filename_xma_sep, device = "png", 
        height = plot_sep_height, ratio = plot_sep_ratio)
