###############################################################################
## FIGURE                                                                    ##
## Xiphophorus maculatus constant-region sashimi plots                       ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")
source("aux/gviz.R")

# Configure input paths
ranges_path <- "../_Data/ranges/xma/xma-new-locus-ranges.tsv"
bam_path_cmd <- "../_Data/bam/xma/xma-new-cmd-vs-PRJNA420092.sorted.bam"
bam_path_cz1 <- "../_Data/bam/xma/xma-new-cz1-vs-PRJNA420092.sorted.bam"
bam_path_cz2 <- "../_Data/bam/xma/xma-new-cz2-vs-PRJNA420092.sorted.bam"

# Configure output
filename <- "xma-new-locus-sashimi"

# Configure locus and cut info
locus_width <- 292642
cut_start <- list(CM = 275000, CD = 275000, CZ1 = 1, CZ2 = 250000)
cut_end <- list(CM = locus_width, CD = locus_width, CZ1 = 8000, CZ2 = 260000)
cut_width <- vapply(names(cut_start), USE.NAMES = TRUE, FUN.VALUE = 0,
                    FUN = function(x) cut_end[[x]] - cut_start[[x]] + 1)
cut_tab <- tibble(feature = names(cut_start), cut_start = as.integer(cut_start),
                  cut_end = as.integer(cut_end), cut_width = cut_width)

#------------------------------------------------------------------------------
# PREPARE CH TABLES
#------------------------------------------------------------------------------

# Get all ranges
ranges_tab <- suppressMessages(read_tsv(ranges_path))

# Subset to CH and extract exon information
ch_tab <- ranges_tab %>% filter(type == "CH") %>%
  mutate(exon = paste0("C", sub("IGH", "", label)),
         feature = sub("-.*", "", exon))

# Split into isoforms by exon lists
isoforms <- list(
  "IGHM-TM" = c("CM-1", "CM-2", "CM-3", "CM-TM1", "CM-TM2"),
  "IGHM-S" = c("CM-1", "CM-2", "CM-3", "CM-4"),
  "IGHD-TM" = c("CM-1", "CD-1", "CD-2A", "CD-3A", "CD-4A",
              "CD-2B", "CD-3B", "CD-4B", "CD-5", "CD-6", "CD-7",
              "CD-TM1", "CD-TM2"),
  "IGHZ1-TM" = c("CZ1-1", "CZ1-2", "CZ1-3", "CZ1-4", "CZ1-TM1", "CZ1-TM2"),
  "IGHZ1-S" = c("CZ1-1", "CZ1-2", "CZ1-3", "CZ1-4", "CZ1-S"),
  "IGHZ2-TM" = c("CZ2-1", "CZ2-2", "CZ2-3", "CZ2-4", "CZ2-TM1", "CZ2-TM2"),
  "IGHZ2-S" = c("CZ2-1", "CZ2-2", "CZ2-3", "CZ2-4", "CZ2-S")
)

ch_tab_isoforms <- lapply(names(isoforms), function(x)
  ch_tab %>% filter(exon %in% isoforms[[x]]) %>% mutate(transcript = x)
) %>% bind_rows %>% inner_join(., cut_tab, by="feature")
ch_tab_cut <- ch_tab_isoforms %>%
  mutate(start = start - cut_start + 1, 
         end = end - cut_start + 1) %>%
  filter(start > 0, end < cut_width)

# Define introns based on isoform
ch_tab_introns <- ch_tab_cut %>%
  group_by(transcript) %>% arrange(start, end) %>%
  mutate(start1 = lag(end) + 1,
         end1 = start - 1,
         start = ifelse(is.na(start1), 1, start1),
         end = ifelse(is.na(end1), cut_width, end1),
         width = end - start + 1,
         strand = "+") %>% ungroup() %>%
  select(-label, -exon, -start1, -end1)

#------------------------------------------------------------------------------
# DEFINE GENOME AXIS TRACKS
#------------------------------------------------------------------------------

# "Chromosome" names for tracks
chr_name <- list(
  CM = "igh_xma_new",
  CD = "igh_xma_new",
  CZ1 = "igh_xma_new",
  CZ2 = "igh_xma_new"
)

# Greate GenomeAxisTrack
xtrack_cm <- GenomeAxisTrack(from=1, to=cut_width[["CM"]], 
                              chromosome=chr_name[["CM"]], 
                          col=colours[["GR"]], fontcolor=colours[["GR"]],
                          size = 8, add35 = FALSE, add53 = FALSE,
                          scale = 2000, labelPos = "beside")

xtrack_cd <- xtrack_cm

xtrack_cz1 <- GenomeAxisTrack(from=1, to=cut_width[["CZ1"]], 
                             chromosome=chr_name[["CZ1"]], 
                          col=colours[["GR"]], fontcolor=colours[["GR"]],
                          size = 8, add35 = FALSE, add53 = FALSE,
                          scale = 2000, labelPos = "beside")

xtrack_cz2 <- GenomeAxisTrack(from=1, to=cut_width[["CZ2"]], 
                             chromosome=chr_name[["CZ2"]], 
                             col=colours[["GR"]], fontcolor=colours[["GR"]],
                             size = 8, add35 = FALSE, add53 = FALSE,
                             scale = 2000, labelPos = "beside")


#------------------------------------------------------------------------------
# DEFINE ALIGNMENT TRACKS
#------------------------------------------------------------------------------

altrack_cut_cmd <- AlignmentsTrack(bam_path_cmd,
                               ispaired = TRUE,
                               type = c("coverage", "sashimi"),
                               chromosome = chr_name[["CM"]],
                               from = 1,
                               to = cut_width[["CM"]],
                               name = "Reads",
                               size = 0.5
)

altrack_cut_cz1 <- AlignmentsTrack(bam_path_cz1,
                                   ispaired = TRUE,
                                   type = c("coverage", "sashimi"),
                                   chromosome = chr_name[["CZ1"]],
                                   from = 1,
                                   to = cut_width[["CZ1"]],
                                   name = "Reads",
                                   size = 0.5
)

altrack_cut_cz2 <- AlignmentsTrack(bam_path_cz2,
                                   ispaired = TRUE,
                                   type = c("coverage", "sashimi"),
                                   chromosome = chr_name[["CZ2"]],
                                   from = 1,
                                   to = cut_width[["CZ2"]],
                                   name = "Reads",
                                   size = 0.5
)

#------------------------------------------------------------------------------
# DEFINE GENE REGION TRACKS
#------------------------------------------------------------------------------

# Define subfigure order
subfigure_labels <- LETTERS
subfigure_order <- seq(length(chr_name))
names(subfigure_order) <- c("CZ1", "CZ2", "CM", "CD")

# Group isoforms by isotype and edit to match alignment range
cm_tab_cut <- ch_tab_cut %>% filter(grepl("IGHM", transcript))
cd_tab_cut <- ch_tab_cut %>% filter(grepl("IGHD", transcript))
cz1_tab_cut <- ch_tab_cut %>% filter(grepl("IGHZ1", transcript)) %>%
  mutate(feature = sub("CZ1", "CZ", feature))
cz2_tab_cut <- ch_tab_cut %>% filter(grepl("IGHZ2", transcript)) %>%
  mutate(feature = sub("CZ2", "CZ", feature))

# Propagate CZ colour to CZ1 and CZ2

grtrack_cut_cm <- GeneRegionTrack(cm_tab_cut,
                                  group = cm_tab_cut$transcript,
                                  feature = cm_tab_cut$feature,
                                  transcript = cm_tab_cut$transcript,
                                  name = subfigure_labels[subfigure_order[["CM"]]],
                                  rotation.title = 0,
                                  col.title = "black",
                                  chromosome = chr_name[["CM"]],
                                  stacking = "squish",
                                  from = 1,
                                  to = cut_width[["CM"]],
                                  size = 10,
                                  background.title="transparent",
                                  transcriptAnnotation = "transcript",
                                  cex.title = 3.5,
                                  fontface.title = 1
)

grtrack_cut_cd <- GeneRegionTrack(cd_tab_cut,
                                  group = cd_tab_cut$transcript,
                                  feature = cd_tab_cut$feature,
                                  transcript = cd_tab_cut$transcript,
                                  name = subfigure_labels[subfigure_order[["CD"]]],
                                  rotation.title = 0,
                                  cex.title = 3.5,
                                  col.title = "black",
                                  fontface.title = 1,
                                  chromosome = chr_name[["CD"]],
                                  stacking = "squish",
                                  from = 1,
                                  to = cut_width[["CD"]],
                                  size = 10,
                                  background.title="transparent",
                                  transcriptAnnotation = "transcript"
)

grtrack_cut_cz1 <- GeneRegionTrack(cz1_tab_cut,
                                  group = cz1_tab_cut$transcript,
                                  feature = cz1_tab_cut$feature,
                                  transcript = cz1_tab_cut$transcript,
                                  name = subfigure_labels[subfigure_order[["CZ1"]]],
                                  rotation.title = 0,
                                  cex.title = 3.5,
                                  col.title = "black",
                                  fontface.title = 1,
                                  chromosome = chr_name[["CZ1"]],
                                  stacking = "squish",
                                  from = 1,
                                  to = cut_width[["CZ1"]],
                                  size = 10,
                                  background.title="transparent",
                                  transcriptAnnotation = "transcript"
)

grtrack_cut_cz2 <- GeneRegionTrack(cz2_tab_cut,
                                   group = cz2_tab_cut$transcript,
                                   feature = cz2_tab_cut$feature,
                                   transcript = cz2_tab_cut$transcript,
                                   name = subfigure_labels[subfigure_order[["CZ2"]]],
                                   rotation.title = 0,
                                   cex.title = 3.5,
                                   col.title = "black",
                                   fontface.title = 1,
                                   chromosome = chr_name[["CZ2"]],
                                   stacking = "squish",
                                   from = 1,
                                   to = cut_width[["CZ2"]],
                                   size = 10,
                                   background.title="transparent",
                                   transcriptAnnotation = "transcript"
)

#------------------------------------------------------------------------------
# DEFINE ISOTYPE CO-ORDINATE LIMITS AND PLOTTING PARAMETERS
#------------------------------------------------------------------------------
title.width <- 1.3

isotype_min <- lapply(tolower(names(subfigure_order)), function(a)
  get(paste0(a, "_tab_cut")) %>% pull(start) %>% min)

isotype_max <- lapply(tolower(names(subfigure_order)), function(a)
  get(paste0(a, "_tab_cut")) %>% pull(end) %>% max)

splice_score <- list(
  CM = 1,#60,
  CD = 10,#30,
  CZ1 = 10,#30
  CZ2 = 10
)
splice_score <- splice_score[match(names(splice_score), names(subfigure_order))]

plot_unit = "cm"

prepare_introns <- function(isotype){
  ch_tab_introns %>% filter(grepl(paste0("IGH", isotype), transcript)) %>%
    mutate(seqnames = chr_name[[paste0("C", isotype)]]) %>%
    select(-transcript) %>% group_by(seqnames, start, end, strand) %>%
    filter(row_number() == 1) %>% 
    as("GRanges")
}

cd_introns <- prepare_introns("D")
cm_introns <- prepare_introns("M")
cz1_introns <- prepare_introns("Z1")
cz2_introns <- prepare_introns("Z2")

#------------------------------------------------------------------------------
# COMBINE AND PLOT TRACKS
#------------------------------------------------------------------------------

# Define flanking regions for plots
flank_start <- list(CM = 700, CD = 1500, CZ1 = 400, CZ2 = 400)
flank_start <- flank_start[match(names(flank_start), names(subfigure_order))]
flank_end <- list(CM = 500, CD = 500, CZ1 = 500, CZ2 = 500)
flank_end <- flank_end[match(names(flank_end), names(subfigure_order))]
filter_tol <- list(CM = 5, CD = 2, CZ1 = 2, CZ2 = 2)
filter_tol <- filter_tol[match(names(filter_tol), names(subfigure_order))]


# Collate and order isotype data
altrack_cut_cm <- altrack_cut_cmd
altrack_cut_cd <- altrack_cut_cmd
grtracks <- lapply(tolower(names(subfigure_order)), function(a)
  get(paste0("grtrack_cut_", a)))
altracks <- lapply(tolower(names(subfigure_order)), function(a)
  get(paste0("altrack_cut_", a)))
xtracks <- lapply(tolower(names(subfigure_order)), function(a)
  get(paste0("xtrack_", a)))
introns <- lapply(tolower(names(subfigure_order)), function(a)
  get(paste0(a, "_introns")))

# Prepare parent viewport
plot_height <- 24.3*1.5
plot_width <- 20.9*1.5
fig_height_ratios <- rep(1, length(grtracks))
map_layout <- split_layout(plot_width, plot_height, 
                           height_ratios = fig_height_ratios)
vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)

# Add subfigures to rows in order
for (n in 1:length(subfigure_order)){
  pushViewport(viewport(layout.pos.col = 1, layout.pos.row = n))
  plotTracks(list(grtracks[[n]], xtracks[[n]], altracks[[n]]),
             from = isotype_min[[n]] - flank_start[[n]],
             to = isotype_max[[n]] + flank_end[[n]],
             sashimiScore = splice_score[[n]], add = TRUE,
             title.width = title.width,
             sashimiFilter = introns[[n]],
             sashimiFilterTolerange = filter_tol[[n]])
  popViewport(1)
}

# TODO: Refine and expand

plt <- grid.grab()
savefig(plot = plt, filename = filename,
        height = plot_height, width = plot_width)