###############################################################################
## FIGURE                                                                    ##
## Xiphophorus maculatus IgH locus map                                       ##
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
locus_path <- "/home/will/Documents/data/igh/other_loci/xma/locus.fa"
vh_ranges_path <- "/home/will/Documents/data/igh/other_loci/xma/vsearch/xma_vh_ranges.tsv"
dh_ranges_path <- "/home/will/Documents/data/igh/other_loci/xma/dsearch/xma_dh_ranges.tsv"
jh_ranges_path <- "/home/will/Documents/data/igh/other_loci/xma/jsearch/xma_jh_ranges.tsv"
ch_ranges_path <- "/home/will/Documents/data/igh/other_loci/xma/csearch/xma_ch_ranges.tsv"

# Configure output
filename <- "xma-locus-map"

#------------------------------------------------------------------------------
# LOCUS IDEOGRAM
#------------------------------------------------------------------------------

# Specify locus and chromosome information
genome_name <- "xma_washu_4.4.2"
chromosome_id <- "16"
chromosome_name <- "16"
chromosome_length <- 25094516
locus_start <- 9384676 + 98855 - 1 # (on chromosome)
locus_end <- 9384676 + 361512 - 1

# Create chromosome table for ideogram
chromosome_info <- tibble(chrom = chromosome_id, chromStart = 0, 
                   chromEnd = chromosome_length,
                   name = chromosome_name, gieStain = "gneg")

# Create ideogram object
itrack <- IdeogramTrack(genome = genome_name, 
                        chromosome = chromosome_id, 
                        name = chromosome_name, 
                        bands = chromosome_info,
                        from = locus_start,
                        end = locus_end)

#------------------------------------------------------------------------------
# STACKED LOCUS MAP
#------------------------------------------------------------------------------

# Set locus and sublocus information
locus_width <- 262658
subloci_labels <- c("IGH")
subloci_seqnames <- c("locus")
subloci_strands <- c("+")

# Get V/D/J/C ranges for annotation
vh_tab <- suppressMessages(read_tsv(vh_ranges_path))
ch_tab <- suppressMessages(read_tsv(ch_ranges_path))
jh_tab <- suppressMessages(read_tsv(jh_ranges_path))
dh_tab <- suppressMessages(read_tsv(dh_ranges_path))

# Split up CH ranges
cd_tab <- ch_tab %>% filter(grepl("IGHD", label))
cm_tab <- ch_tab %>% filter(grepl("IGHM", label))
cz_tab <- ch_tab %>% filter(grepl("IGHZ", label))

# Get sublocus boundaries
all_tab <- bind_rows(vh_tab, jh_tab, ch_tab, dh_tab)
subloci <- tibble()
for (n in seq(length(subloci_labels))){
  l <- subloci_labels[n]
  t <- filter(all_tab, grepl(l, label)) %>%
    summarise(label = l, start = min(start), end = max(end), 
              width = end - start + 1, 
              seqnames = subloci_seqnames[n],
              strand = subloci_strands[n])
  subloci <- bind_rows(subloci, t)
}

# Greate GenomeAxisTrack
xtrack1 <- GenomeAxisTrack(from=1, to=locus_width, chromosome="locus", 
                           col=colours[["HL1"]], fontcolor=colours[["HL1"]])
xtrack2 <- GenomeAxisTrack(from=1, to=locus_width, chromosome="locus", 
                           col=colours[["HL2"]], fontcolor=colours[["HL2"]],
                           size = 1.1, add35 = FALSE, add53 = FALSE)

# Create AnnotationTrack objects
sublocus_track <- AnnotationTrack(subloci, feature = "Sublocus", name = "Sublocus", 
                                  group = subloci$label, cex = 1.4, shape = "arrow",
                                  featureAnnotation = "group", size=1.5, chromosome = chromosome_id)
vh_track <- AnnotationTrack(vh_tab, feature = "VH", name = "VH", chromosome = chromosome_id,
                            background.title = colours[["VH"]],
                            background.panel = paste0(colours[["VH"]], "20"))
jh_track <- AnnotationTrack(jh_tab, feature = "JH", name = "JH", chromosome = chromosome_id,
                            background.title = colours[["JH"]], 
                            background.panel = paste0(colours[["JH"]], "20"))
dh_track <- AnnotationTrack(dh_tab, feature = "DH", name = "DH", chromosome = chromosome_id,
                            background.title = colours[["DH"]], 
                            background.panel = paste0(colours[["DH"]], "20"))
cz_track <- AnnotationTrack(cz_tab, feature = "CZ", name = "IGHZ", chromosome = chromosome_id,
                            background.title = colours[["CZ"]], 
                            background.panel = paste0(colours[["CZ"]], "20"))
cm_track <- AnnotationTrack(cm_tab, feature = "CM", name = "IGHM", chromosome = chromosome_id,
                            background.title = colours[["CM"]], 
                            background.panel = paste0(colours[["CM"]], "20"))
cd_track <- AnnotationTrack(cd_tab, feature = "CD", name = "IGHD", chromosome = chromosome_id,
                            background.title = colours[["CD"]], 
                            background.panel = paste0(colours[["CD"]], "20"))

#------------------------------------------------------------------------------
# ZOOMED CM/D-REGION MAPS
#------------------------------------------------------------------------------

jgap_max <- 10000

# Re-annotate C-range table and split into subloci
cmd_tab <- ch_tab %>%
  filter(grepl("IGH[MD]", label)) %>%
  mutate(exon = paste0("C", sub("IGH", "", label)),
         feature = sub("-.*", "", exon))

jmd_tab <- jh_tab %>%
  mutate(feature = "JH",
         exon = paste0("JH-", sub("IGHJ", "", label))) %>%
  filter(start >= min(cmd_tab$start) - jgap_max,
         start <= max(cmd_tab$end))


chm_tm_tab <- cmd_tab %>% 
  filter(exon %in% c("CM-1", "CM-2", "CM-3", "CM-TM1", "CM-TM2")) %>%
  mutate(transcript = "IGHM-TM")

chm_s_tab <- cmd_tab %>% 
  filter(exon %in% c("CM-1", "CM-2", "CM-3", "CM-4")) %>%
  mutate(transcript = "IGHM-S")

chd_tm_tab <- cmd_tab %>%
  filter(exon == "CM-1" | grepl("CD", exon)) %>%
  mutate(transcript = "IGHD-TM")

cmd_tab_summary <- bind_rows(cmd_tab, jmd_tab)
cmd_tab_grouped <- bind_rows(chm_tm_tab, chm_s_tab, chd_tm_tab)

# Make overall and transcript-specific tracks

cmd_summary_track <-
  AnnotationTrack(cmd_tab_summary,
                  feature = cmd_tab_summary$feature,
                  group = cmd_tab_summary$exon,
                  chromosome = chromosome_id,
                  groupAnnotation = "group",
                  name = "Exons",
                  stackHeight=0.4,
                  size = 3) # Give enough height to manually add labels later

cmd_grouped_track <-
  GeneRegionTrack(cmd_tab_grouped,
                  group = cmd_tab_grouped$transcript,
                  feature = cmd_tab_grouped$feature,
                  chromosome= chromosome_id,
                  transcript = cmd_tab_grouped$transcript,
                  transcriptAnnotation = "transcript",
                  cex.group = 1,
                  name = "Isoforms")

#------------------------------------------------------------------------------
# ZOOMED CZ-REGION MAPS
#------------------------------------------------------------------------------

cz_tab <- ch_tab %>%
  filter(grepl("IGHZ", label)) %>%
  mutate(exon = paste0("C", sub("IGH", "", label)),
         feature = sub("-.*", "", exon))

jz_tab <- jh_tab %>%
  mutate(feature = "JH",
         exon = paste0("JH-", sub("IGHJ", "", label))) %>%
  filter(start >= min(cz_tab$start) - jgap_max,
         start <= max(cz_tab$end))

dz_tab <- dh_tab %>%
  mutate(feature = "DH",
         exon = paste0("DH-", sub("IGHD", "", label))) %>%
  filter(start >= min(cz_tab$start) - jgap_max,
         start <= max(cz_tab$end))

vz_tab <- vh_tab %>%
  mutate(feature = "VH",
         exon = paste0("VH-", sub("IGHV", "", label))) %>%
  filter(start >= min(cz_tab$start) - jgap_max,
         start <= max(cz_tab$end))

chz_tm_tab <- cz_tab %>% 
  filter(exon %in% c("CZ-1", "CZ-2", "CZ-3", "CZ-4", "CZ-TM1", "CZ-TM2")) %>%
  mutate(transcript = "IGHZ-TM")

chz_s_tab <- cz_tab %>% 
  filter(exon %in% c("CZ-1", "CZ-2", "CZ-3", "CZ-4", "CZ-S")) %>%
  mutate(transcript = "IGHZ-S")

cz_tab_summary <- bind_rows(cz_tab, jz_tab, dz_tab, vz_tab)
cz_tab_grouped <- bind_rows(chz_tm_tab, chz_s_tab)

# Make overall and transcript-specific tracks
cz_summary_track <-
  AnnotationTrack(cz_tab_summary,
                  feature = cz_tab_summary$feature,
                  group = cz_tab_summary$exon,
                  chromosome= chromosome_id,
                  groupAnnotation = "group",
                  name = "Exons",
                  stackHeight=0.4,
                  size = 3)

cz_grouped_track <-
  GeneRegionTrack(cz_tab_grouped,
                  group = cz_tab_grouped$transcript,
                  feature = cz_tab_grouped$feature,
                  chromosome= chromosome_id,
                  transcript = cz_tab_grouped$transcript,
                  transcriptAnnotation = "transcript",
                  name = "Isoforms",
                  size = 1.5)

#------------------------------------------------------------------------------
# COMBINE AND PLOT TRACKS
#------------------------------------------------------------------------------

# Add C-region highlighting to stacked map
htrack <- 
  HighlightTrack(
    trackList = list(sublocus_track, vh_track, dh_track, jh_track, cm_track, cd_track, cz_track),
    start = c(min(cz_tab_summary$start, cz_tab_summary$end)-2000,
              min(cmd_tab_summary$start, cmd_tab_summary$end)-2000),
    end = c(max(cz_tab_summary$start, cz_tab_summary$end)+2000,
            max(cmd_tab_summary$start, cmd_tab_summary$end)+2000),
    chromosome = chromosome_id,
    from = 1,
    to = locus_width
    )

# Combine C-region tracks
ctrack_md <-
  HighlightTrack(
    trackList = list(xtrack2, cmd_summary_track, cmd_grouped_track),
    chromosome = "chr6",
    from = min(cmd_tab_summary$start, cmd_tab_summary$end)-2000,
    to = max(cmd_tab_summary$start, cmd_tab_summary$end)+2000
  )
ctrack_z <-
  HighlightTrack(
    trackList = list(xtrack2, cz_summary_track, cz_grouped_track),
    chromosome = "chr6",
    from = min(cz_tab_summary$start, cz_tab_summary$end)-2000,
    to = max(cz_tab_summary$start, cz_tab_summary$end)+2000
  )


# Stack plots together
width_md <- max(cmd_tab_summary$start, cmd_tab_summary$end)+2000 - (min(cmd_tab_summary$start, cmd_tab_summary$end)-2000)
width_z <- max(cz_tab_summary$start, cz_tab_summary$end)+2000 - (min(cz_tab_summary$start, cz_tab_summary$end)-2000)
width_ratio <- width_md/width_z
width_partition <- c(1, width_ratio)
# Stack plots together
plot_unit = "cm"
plot_height <- 24.3*1.5
plot_width <- 20.9*1.5
fig_height_ratios <- c(1,6,4)
fig_heights <- sapply(fig_height_ratios, 
                      function(x) x/sum(fig_height_ratios) * plot_height)
fig_widths <- sapply(width_partition,
                     function(x) x/sum(width_partition) * plot_width)
map_layout <- grid.layout(
  ncol = 2,
  nrow = 3,
  heights = unit(fig_heights, plot_unit),
  widths = unit(fig_widths, plot_unit)
)
vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
# Add ideogram to first row
pushViewport(viewport(layout.pos.col = c(1,2), layout.pos.row = 1))
plotTracks(itrack, from=locus_start, to=locus_end, add=TRUE)
popViewport(1)
# Add stacked map to second row
pushViewport(viewport(layout.pos.col = c(1,2), layout.pos.row = 2))
plotTracks(list(xtrack1, htrack), add=TRUE)
popViewport(1)
# Add C-region maps to third row
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3))
plotTracks(ctrack_md, add=TRUE)
popViewport(1)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3))
plotTracks(ctrack_z, add=TRUE)
popViewport(1)

plt <- grid.grab()
savefig(plot = plt, filename = filename, device = "svg", 
        height = plot_height, width = plot_width)