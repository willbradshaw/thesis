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
ranges_path <- "../_Data/ranges/xma/xma-locus-ranges.tsv"
bam_path_cmd <- "../_Data/bam/xma/xma-cmd-vs-PRJNA420092.sorted.bam"
bam_path_cz <- "../_Data/bam/xma/xma-cz-vs-PRJNA420092.sorted.bam"

# Configure output
filename <- "xma-locus-sashimi"

# Configure locus and cut info
locus_width <- 262658
cut_start <- list(CM = 245000, CD = 245000, CZ = 1)
cut_end <- list(CM = locus_width, CD = locus_width, CZ = 8000)
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
  "IGHZ-TM" = c("CZ-1", "CZ-2", "CZ-3", "CZ-4", "CZ-TM1", "CZ-TM2"),
  "IGHZ-S" = c("CZ-1", "CZ-2", "CZ-3", "CZ-4", "CZ-S")
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
  CM = "igh_xma",
  CD = "igh_xma",
  CZ = "igh_xma_ighz"
)

# Greate GenomeAxisTrack
xtrack_cmd <- GenomeAxisTrack(from=1, to=cut_width[["CM"]], 
                              chromosome=chr_name[["CM"]], 
                          col=colours[["GR"]], fontcolor=colours[["GR"]],
                          size = 8, add35 = FALSE, add53 = FALSE,
                          scale = 2000, labelPos = "beside")

xtrack_cz <- GenomeAxisTrack(from=1, to=cut_width[["CZ"]], 
                             chromosome=chr_name[["CZ"]], 
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

altrack_cut_cz <- AlignmentsTrack(bam_path_cz,
                                   ispaired = TRUE,
                                   type = c("coverage", "sashimi"),
                                   chromosome = chr_name[["CZ"]],
                                   from = 1,
                                   to = cut_width[["CZ"]],
                                   name = "Reads",
                                   size = 0.5
)

#------------------------------------------------------------------------------
# DEFINE GENE REGION TRACKS
#------------------------------------------------------------------------------

# Group isoforms by isotype and edit to match alignment range
cm_tab_cut <- ch_tab_cut %>% filter(grepl("IGHM", transcript))
cd_tab_cut <- ch_tab_cut %>% filter(grepl("IGHD", transcript))
cz_tab_cut <- ch_tab_cut %>% filter(grepl("IGHZ", transcript))

grtrack_cut_cm <- GeneRegionTrack(cm_tab_cut,
                                  group = cm_tab_cut$transcript,
                                  feature = cm_tab_cut$feature,
                                  transcript = cm_tab_cut$transcript,
                                  name = "A",
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
                                  name = "B",
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

grtrack_cut_cz <- GeneRegionTrack(cz_tab_cut,
                                  group = cz_tab_cut$transcript,
                                  feature = cz_tab_cut$feature,
                                  transcript = cz_tab_cut$transcript,
                                  name = "C",
                                  rotation.title = 0,
                                  cex.title = 3.5,
                                  col.title = "black",
                                  fontface.title = 1,
                                  chromosome = chr_name[["CZ"]],
                                  stacking = "squish",
                                  from = 1,
                                  to = cut_width[["CZ"]],
                                  size = 10,
                                  background.title="transparent",
                                  transcriptAnnotation = "transcript"
)

#------------------------------------------------------------------------------
# DEFINE ISOTYPE CO-ORDINATE LIMITS AND PLOTTING PARAMETERS
#------------------------------------------------------------------------------
title.width <- 1.3

isotype_min <- list(
  CM = cm_tab_cut %>% pull(start) %>% min,
  CD = cd_tab_cut %>% pull(start) %>% min,
  CZ = cz_tab_cut %>% pull(start) %>% min
)

isotype_max <- list(
  CM = cm_tab_cut %>% pull(end) %>% max,
  CD = cd_tab_cut %>% pull(end) %>% max,
  CZ = cz_tab_cut %>% pull(end) %>% max
)

splice_score <- list(
  CM = 1,#60,
  CD = 10,#30,
  CZ = 10#30
)

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
cz_introns <- prepare_introns("Z")

#------------------------------------------------------------------------------
# COMBINE AND PLOT TRACKS
#------------------------------------------------------------------------------

# Prepare parent viewport
plot_height <- 24.3*1.5
plot_width <- 20.9*1.5
fig_height_ratios <- c(1,1,1)
map_layout <- split_layout(plot_width, plot_height, 
                           height_ratios = fig_height_ratios)
vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
# Add CM to first row
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
plotTracks(list(grtrack_cut_cm, xtrack, altrack_cut_cmd), 
           from = isotype_min[["CM"]]-700, to = isotype_max[["CM"]]+500, 
           sashimiScore = splice_score[["CM"]], add = TRUE,
           title.width = title.width,
           sashimiFilter = cm_introns,
           sashimiFilterTolerance = 5)
popViewport(1)
# Add CM1/CD to second row
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
plotTracks(list(grtrack_cut_cd, xtrack, altrack_cut_cmd), 
           from = isotype_min[["CM"]]-1500, to = isotype_max[["CD"]]+500, 
           sashimiScore = splice_score[["CD"]], add = TRUE,
           title.width = title.width,
           sashimiFilter = cd_introns,
           sashimiFilterTolerance = 2)
# plotTracks(list(grtrack_cut_cd, xtrack, altrack_cut_cmd), 
#            from = cm_min-1800, to = cd_max+500, 
#            sashimiFilter = cd_introns, add = TRUE,
#            title.width = title.width)
popViewport(1)
# Add CZ to third row
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3))
plotTracks(list(grtrack_cut_cz, xtrack, altrack_cut_cz), 
           from = isotype_min[["CZ"]]-400, to = isotype_max[["CZ"]]+500, 
           sashimiScore = splice_score[["CZ"]], add = TRUE,
           title.width = title.width,
           sashimiFilter = cz_introns,
           sashimiFilterTolerance = 2) # TODO: Modify to include V/D/J?
popViewport(1)


# TODO: Refine and expand

plt <- grid.grab()
savefig(plot = plt, filename = filename, device = "svg", 
        height = plot_height, width = plot_width)
savefig(plot = plt, filename = filename, device = "png", 
        height = plot_height, width = plot_width)

