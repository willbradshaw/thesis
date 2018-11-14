###############################################################################
## FIGURE                                                                    ##
## Nothobranchius furzeri RSS composition plots                              ##
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
vh_rss_path_nfu <- "../_Data/rss/nfu/nfu_vh_rss.fasta"
dh_rss_path_nfu <- "../_Data/rss/nfu/nfu_dh_rss.fasta"
jh_rss_path_nfu <- "../_Data/rss/nfu/nfu_jh_rss.fasta"

# Configure output
filename_nfu_all <- "nfu-rss-seqlogo-all"
filename_nfu_sep <- "nfu-rss-seqlogo-sep"

#------------------------------------------------------------------------------
# IMPORT RSS SEQUENCES AND EXTRACT HEPTAMER/NONAMER/SPACER SEQUENCES
#------------------------------------------------------------------------------

# Get D/J sequences
vh_rss_nfu <- readDNAStringSet(vh_rss_path_nfu)
dh_rss_nfu <- readDNAStringSet(dh_rss_path_nfu)
jh_rss_nfu <- readDNAStringSet(jh_rss_path_nfu)

# Group by length
rss_nfu <- list(
  vh = vh_rss_nfu,
  dh = dh_rss_nfu,
  jh = jh_rss_nfu,
  long = c(vh_rss_nfu, jh_rss_nfu),
  short = dh_rss_nfu,
  all = c(vh_rss_nfu, dh_rss_nfu, jh_rss_nfu)
)

# Extract subsequences

heptamers_nfu <- lapply(rss_nfu, function(seqs) subseq(seqs, 1, 7))
nonamers_nfu <- lapply(rss_nfu, function(seqs) subseq(seqs, -9, -1))
spacers_nfu <- lapply(rss_nfu, function(seqs) subseq(seqs, 8, -10))

#------------------------------------------------------------------------------
# COMPUTE PWMS AND LENGTH DISTRIBUTIONS
#------------------------------------------------------------------------------

# Character vectors
charHeptamer_nfu <- lapply(heptamers_nfu, function(seqs) as.character(seqs))
charNonamer_nfu <- lapply(nonamers_nfu, function(seqs) as.character(seqs))

# Spacer length tables
getLengthTable <- function(stringSet, type, length){
  tibble(name = names(stringSet), width = width(stringSet),
                seq = as.character(stringSet), type=type, length=length)
}
lengths <- list(vh="long", dh="short", jh="long", long="long", short="short", all="all")
spacerTable_nfu <- bind_rows(
  getLengthTable(spacers_nfu$vh, "VH", "long"),
  getLengthTable(spacers_nfu$dh, "DH", "short"),
  getLengthTable(spacers_nfu$jh, "JH", "long")
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
ggHeptamer_nfu <- lapply(charHeptamer_nfu, function(seqs) ggplot() + gglogo(seqs) +
                       xlab("Position") + theme_logo)
ggNonamer_nfu <- lapply(charNonamer_nfu, function(seqs) ggplot() + gglogo(seqs) +
                      xlab("Position") + theme_logo)


# Length distributions
tick_colours = rep(c("black", "red", "black"),2)

ggSpacer_nfu <- ggplot(spacerTable_nfu) + 
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
annot_coords_nfu <- list(
  VH = list(x = 19.8, y = 14),
  DH = list(x = 14.2, y = 23),
  JH = list(x = 25, y = 12)
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

spacer_all_annotated_nfu <- annotate_spacer_plot(ggSpacer_nfu, annot_coords_nfu)

plt_nfu_all <- plot_grid(ggHeptamer_nfu$all + ggtitle("Heptamer composition"),
             spacer_all_annotated_nfu + ggtitle("Spacer length"),
             ggNonamer_nfu$all + ggtitle("Nonamer composition"),
             ncol = 3, labels="AUTO", label_fontfamily = titlefont,
             label_fontface = "plain",
             label_size = fontsize_base * fontscale_label)

#------------------------------------------------------------------------------
# SPLIT-BY-SEGMENT-TYPE SEQLOGO PLOT
#------------------------------------------------------------------------------

plt_nfu_sep <- plot_grid(ggHeptamer_nfu$vh + ggtitle("VH heptamers"),
                         ggNonamer_nfu$vh + ggtitle("VH nonamers"),
                         ggHeptamer_nfu$dh + ggtitle("DH heptamers"),
                         ggNonamer_nfu$dh + ggtitle("DH nonamers"),
                         ggHeptamer_nfu$jh + ggtitle("JH heptamers"),
                         ggNonamer_nfu$jh + ggtitle("JH nonamers"),
                         ncol = 2, labels="AUTO", label_fontfamily = titlefont,
                         label_fontface = "plain",
                         label_size = fontsize_base * fontscale_label)

#------------------------------------------------------------------------------
# SAVE PLOTS
#------------------------------------------------------------------------------

plot_all_height <- 10
plot_all_ratio <- 3 # Width vs height

savefig(plot = plt_nfu_all, filename = filename_nfu_all, device = "svg", 
        height = plot_all_height, ratio = plot_all_ratio)
savefig(plot = plt_nfu_all, filename = filename_nfu_all, device = "png", 
        height = plot_all_height, ratio = plot_all_ratio)

plot_sep_height <- 30
plot_sep_ratio <- 2/3 # Width vs height

savefig(plot = plt_nfu_sep, filename = filename_nfu_sep, device = "svg", 
        height = plot_sep_height, ratio = plot_sep_ratio)
savefig(plot = plt_nfu_sep, filename = filename_nfu_sep, device = "png", 
        height = plot_sep_height, ratio = plot_sep_ratio)
