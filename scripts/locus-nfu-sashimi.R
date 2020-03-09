###############################################################################
## FIGURE                                                                    ##
## Nothobranchius furzeri constant-region sashimi plots                      ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# PREPARE CM TABLES
#------------------------------------------------------------------------------

# Configure Gviz with locus info
chromosome_id <- "locus"
locus_width <- 306408
subloci_labels <- c("IGH1", "IGH2")
subloci_seqnames <- c("locus", "locus")
subloci_strands <- c("+", "-")

# Re-annotate C-range table and split into subloci
ch_tab <- suppressMessages(read_tsv(inpath_ch)) %>%
  mutate(exon = paste0("C", sub("IGH\\d", "", label)),
         feature = sub("-.*", "", exon))
ch_tab_subloci <- lapply(subloci_labels, function(l) ch_tab %>% filter(grepl(l, label)))
names(ch_tab_subloci) <- subloci_labels

chm_tm_tab <- lapply(subloci_labels, function(l)
  ch_tab_subloci[[l]] %>% filter(exon %in% c("CM-1", "CM-2", "CM-TM1", "CM-TM2")) %>%
    mutate(transcript = "IGHM-TM"))
names(chm_tm_tab) <- subloci_labels

chm_s_tab <- lapply(subloci_labels, function(l)
  ch_tab_subloci[[l]] %>% filter(exon %in% c("CM-1", "CM-2", "CM-3", "CM-4")) %>%
    mutate(transcript = "IGHM-S"))
names(chm_s_tab) <- subloci_labels

chd_tm_tab <- lapply(subloci_labels, function(l)
  ch_tab_subloci[[l]] %>% filter(exon == "CM-1" | grepl("CD", exon)) %>%
    mutate(transcript = "IGHD-TM"))
names(chd_tm_tab) <- subloci_labels

ch_tab_subloci_grouped <- lapply(subloci_labels, function(l)
  bind_rows(chm_tm_tab[[l]], chm_s_tab[[l]], chd_tm_tab[[l]]))
names(ch_tab_subloci_grouped) <- subloci_labels

# Edit to match alignment range and split by isotype
cut_start <- 128001
cut_end <- 149000
cut_width <- cut_end - cut_start + 1

ch_tab_subloci_cut <- ch_tab_subloci_grouped[["IGH1"]] %>%
  mutate(start = start - cut_start + 1,
         end = end - cut_start + 1) %>%
  filter(start > 0, end < cut_width)

cm_tab_subloci_cut <- ch_tab_subloci_cut %>% 
  filter(grepl("IGHM", transcript))
cd_tab_subloci_cut <- ch_tab_subloci_cut %>% 
  filter(grepl("IGHD", transcript))

# Define CD-introns (complement of CD-exons)
cd_introns_tab <- cd_tab_subloci_cut %>% 
  arrange(start, end) %>%
  mutate(start1 = lag(end) + 1,
         end1 = start - 1,
         start = ifelse(is.na(start1), 1, start1),
         end = ifelse(is.na(end1), cut_width, end1),
         width = end - start + 1,
         strand = "+") %>%
  select(-label, -exon, -feature, -start1, -end1)

cd_introns <- cd_introns_tab %>% as("GRanges")

# Greate GenomeAxisTrack
xtrack <- GenomeAxisTrack(from=1, to=locus_width, chromosome=chromosome_id,
                           col=colours[["GR"]], fontcolor=colours[["GR"]],
                           size = 8, add35 = FALSE, add53 = FALSE,
                           scale = 2000, labelPos = "beside")


#------------------------------------------------------------------------------
# CM ALIGNMENT TRACK
#------------------------------------------------------------------------------
title.width <- 1.3

altrack_cut <- AlignmentsTrack(inpath_bam,
                               ispaired = TRUE,
                               type = c("coverage", "sashimi"),
                               chromosome = chromosome_id,
                               from = 1,
                               to = cut_width,
                               name = "Reads",
                               size = 0.5,
                               )


grtrack_cut_cm <- GeneRegionTrack(cm_tab_subloci_cut,
                               group = cm_tab_subloci_cut$transcript,
                               feature = cm_tab_subloci_cut$feature,
                               transcript = cm_tab_subloci_cut$transcript,
                               name = "A",
                               rotation.title = 0,
                               col.title = "black",
                               chromosome = chromosome_id,,
                               stacking = "squish",
                               from = 1,
                               to = cut_width,
                               size = 10,
                               #stackHeight = 1,
                               background.title="transparent",
                               transcriptAnnotation = "transcript",
                               cex.title = 3.5,
                               fontface.title = 1,
)

cm_min <- ch_tab_subloci_cut %>% filter(feature == "CM") %>% pull(start) %>% min
cm_max <- ch_tab_subloci_cut %>% filter(feature == "CM") %>% pull(end) %>% max
cm_score <- 60

#------------------------------------------------------------------------------
# CD ALIGNMENT TRACK
#------------------------------------------------------------------------------

grtrack_cut_cd <- GeneRegionTrack(cd_tab_subloci_cut,
                                  group = cd_tab_subloci_cut$transcript,
                                  feature = cd_tab_subloci_cut$feature,
                                  transcript = cd_tab_subloci_cut$transcript,
                                  name = "B",
                                  rotation.title = 0,
                                  cex.title = 3.5,
                                  col.title = "black",
                                  fontface.title = 1,
                                  chromosome = chromosome_id,
                                  stacking = "squish",
                                  from = 1,
                                  to = cut_width,
                                  size = 10,
                                  background.title="transparent",
                                  transcriptAnnotation = "transcript",
                                  )

cd_min <- ch_tab_subloci_cut %>% filter(feature == "CD") %>% pull(start) %>% min
cd_max <- ch_tab_subloci_cut %>% filter(feature == "CD") %>% pull(end) %>% max
cd_score <- 30

#------------------------------------------------------------------------------
# COMBINE AND PLOT TRACKS
#------------------------------------------------------------------------------

# Stack plots together
plot_unit = "cm"
plot_height <- 24.3*1.2
plot_width <- 20.9*1.5
fig_height_ratios <- c(1,1)
fig_heights <- sapply(fig_height_ratios, 
                      function(x) x/sum(fig_height_ratios) * plot_height)
fig_widths <- plot_width
map_layout <- grid.layout(
  ncol = length(fig_widths),
  nrow = length(fig_heights),
  heights = unit(fig_heights, plot_unit),
  widths = unit(fig_widths, plot_unit)
)

vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
# Add CM to first row
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
plotTracks(list(grtrack_cut_cm, xtrack, altrack_cut), 
           from = cm_min-800, to = cm_max+500, 
           sashimiScore = cm_score, add = TRUE,
           title.width = title.width)
popViewport(1)
# # Add CD to second row
# pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
# plotTracks(list(grtrack_cut_cd, altrack_cut), 
#            from = cd_min-500, to = cd_max+500, 
#            sashimiScore = cd_score, add = TRUE,
#            title.width = title.width)
# popViewport(1)
# Add CM1/CD to second row
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
plotTracks(list(grtrack_cut_cd, xtrack, altrack_cut), 
           from = cm_min-1800, to = cd_max+500, 
           sashimiFilter = cd_introns, add = TRUE,
           title.width = title.width)
popViewport(1)

# TODO: Refine and expand

plt <- grid.grab()

ggsave(plot = plt, filename = outpath, device = "svg", height=plot_height,
       width = plot_width, units = "cm")
