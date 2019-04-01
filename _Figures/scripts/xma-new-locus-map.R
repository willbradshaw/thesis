###############################################################################
## FIGURE                                                                    ##
## Xiphophorus maculatus IgH locus map (new locus)                           ##
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

# Install Gviz version with tick-decreasing option
gvizV1 <- packageVersion("Gviz")
install.packages("http://bioconductor.org/packages/devel/bioc/src/contrib/Gviz_1.27.6.tar.gz", repos=NULL)
library(Gviz)
if (packageVersion("Gviz") != "1.27.6") stop("Wrong Gviz version loaded.")

# Configure input paths
locus_path <- "../_Data/locus/complete/xma_new.fasta"
vh_ranges_path <- "../_Data/ranges/xma/xma-new-vh-ranges.tsv"
dh_ranges_path <- "../_Data/ranges/xma/xma-new-dh-ranges.tsv"
jh_ranges_path <- "../_Data/ranges/xma/xma-new-jh-ranges.tsv"
ch_ranges_path <- "../_Data/ranges/xma/xma-new-ch-ranges.tsv"

# Configure output
filename <- "xma-new-locus-map"

#------------------------------------------------------------------------------
# LOCUS IDEOGRAM
#------------------------------------------------------------------------------

# Specify locus and chromosome information
genome_name <- "X_maculatus-5.0-male"
chromosome_id <- "16"
chromosome_name <- "16"
chromosome_length <- 25766145
locus_start <- 15758158 # (on chromosome)
locus_end <- 16050799 # NB: Locus is in antisense on this (new) genome assembly

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
locus_width <- 292642
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
                           col=colours[["HL1"]], fontcolor=colours[["HL1"]],
                           add35 = FALSE, add53 = FALSE)
xtrack2 <- GenomeAxisTrack(from=1, to=locus_width, chromosome="locus", 
                           col=colours[["HL2"]], fontcolor=colours[["HL2"]],
                           size = 1.1, add35 = FALSE, add53 = FALSE,
                           distFromAxis=4)

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

jgap_max <- 5000

# Re-annotate C-range table and split into subloci
cmd_tab <- ch_tab %>%
  filter(grepl("IGH[MD]", label)) %>%
  mutate(exon1 = paste0("C", sub("IGH", "", label)),
         feature = sub("-.*", "", exon1))

jmd_tab <- jh_tab %>%
  mutate(feature = "JH",
         exon1 = paste0("JH-", sub("IGHJ", "", label))) %>%
  filter(start >= min(cmd_tab$start) - jgap_max,
         start <= max(cmd_tab$end))


chm_tm_tab <- cmd_tab %>% 
  filter(exon1 %in% c("CM-1", "CM-2", "CM-3", "CM-TM1", "CM-TM2")) %>%
  mutate(transcript = "IGHM-TM")

chm_s_tab <- cmd_tab %>% 
  filter(exon1 %in% c("CM-1", "CM-2", "CM-3", "CM-4")) %>%
  mutate(transcript = "IGHM-S")

chd_tm_tab <- cmd_tab %>%
  filter(exon1 == "CM-1" | grepl("CD", exon1)) %>%
  mutate(transcript = "IGHD-TM")

cmd_tab_summary <- bind_rows(cmd_tab, jmd_tab)
cmd_tab_grouped <- bind_rows(chm_tm_tab, chm_s_tab, chd_tm_tab)

# Make overall and transcript-specific tracks

cmd_summary_track <-
  AnnotationTrack(cmd_tab_summary,
                  feature = cmd_tab_summary$feature,
                  group = cmd_tab_summary$exon1,
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

# IGHZ1

cz1_tab <- ch_tab %>%
  filter(grepl("IGHZ1", label)) %>%
  mutate(exon1 = paste0("C", sub("IGHZ1", "Z", label)),
         feature = sub("-.*", "", exon1))

jz1_tab <- jh_tab %>%
  mutate(feature = "JH",
         exon1 = paste0("JH-", sub("IGHJ", "", label))) %>%
  filter(start >= min(cz1_tab$start) - jgap_max,
         start <= max(cz1_tab$end))

dz1_tab <- dh_tab %>%
  mutate(feature = "DH",
         exon1 = paste0("DH-", sub("IGHD", "", label))) %>%
  filter(start >= min(cz1_tab$start) - jgap_max,
         start <= max(cz1_tab$end))

vz1_tab <- vh_tab %>%
  mutate(feature = "VH",
         exon1 = paste0("VH-", sub("IGHV", "", label))) %>%
  filter(start >= min(cz1_tab$start) - jgap_max,
         start <= max(cz1_tab$end))

cz1_tm_tab <- cz1_tab %>% 
  filter(exon1 %in% c("CZ-1", "CZ-2", "CZ-3", "CZ-4", "CZ-TM1", "CZ-TM2")) %>%
  mutate(transcript = "IGHZ1-TM")

cz1_s_tab <- cz1_tab %>% 
  filter(exon1 %in% c("CZ-1", "CZ-2", "CZ-3", "CZ-4", "CZ-S")) %>%
  mutate(transcript = "IGHZ1-S")

cz1_tab_summary <- bind_rows(cz1_tab, jz1_tab, dz1_tab, vz1_tab)
cz1_tab_grouped <- bind_rows(cz1_tm_tab, cz1_s_tab)

# Make overall and transcript-specific tracks
cz1_summary_track <-
  AnnotationTrack(cz1_tab_summary,
                  feature = cz1_tab_summary$feature,
                  group = cz1_tab_summary$exon1,
                  chromosome= chromosome_id,
                  groupAnnotation = "group",
                  name = "Exons",
                  stackHeight=0.4,
                  size = 3)

cz1_grouped_track <-
  GeneRegionTrack(cz1_tab_grouped,
                  group = cz1_tab_grouped$transcript,
                  feature = cz1_tab_grouped$feature,
                  chromosome= chromosome_id,
                  transcript = cz1_tab_grouped$transcript,
                  transcriptAnnotation = "transcript",
                  name = "Isoforms",
                  size = 1.5)

# IGHZ2

cz2_tab <- ch_tab %>%
  filter(grepl("IGHZ2", label)) %>%
  mutate(exon1 = paste0("C", sub("IGHZ2", "Z", label)),
         feature = sub("-.*", "", exon1))

jz2_tab <- jh_tab %>%
  mutate(feature = "JH",
         exon1 = paste0("JH-", sub("IGHJ", "", label))) %>%
  filter(start >= min(cz2_tab$start) - jgap_max,
         start <= max(cz2_tab$end))

dz2_tab <- dh_tab %>%
  mutate(feature = "DH",
         exon1 = paste0("DH-", sub("IGHD", "", label))) %>%
  filter(start >= min(cz2_tab$start) - jgap_max,
         start <= max(cz2_tab$end))

cz2_tm_tab <- cz2_tab %>% 
  filter(exon1 %in% c("CZ-1", "CZ-2", "CZ-3", "CZ-4", "CZ-TM1", "CZ-TM2")) %>%
  mutate(transcript = "IGHZ2-TM")

# cz2_s_tab <- cz2_tab %>% 
#   filter(exon1 %in% c("CZ-1", "CZ-2", "CZ-3", "CZ-4", "CZ-S")) %>%
#   mutate(transcript = "IGHZ2-S")

cz2_tab_summary <- bind_rows(cz2_tab, jz2_tab, dz2_tab)
cz2_tab_grouped <- bind_rows(cz2_tm_tab)#, cz2_s_tab)

# Make overall and transcript-specific tracks
cz2_summary_track <-
  AnnotationTrack(cz2_tab_summary,
                  feature = cz2_tab_summary$feature,
                  group = cz2_tab_summary$exon1,
                  chromosome= chromosome_id,
                  groupAnnotation = "group",
                  name = "Exons",
                  stackHeight=0.4,
                  size = 3)

cz2_grouped_track <-
  GeneRegionTrack(cz2_tab_grouped,
                  group = cz2_tab_grouped$transcript,
                  feature = cz2_tab_grouped$feature,
                  chromosome= chromosome_id,
                  transcript = cz2_tab_grouped$transcript,
                  transcriptAnnotation = "transcript",
                  name = "Isoforms",
                  size = 3)

#------------------------------------------------------------------------------
# COMBINE AND PLOT TRACKS
#------------------------------------------------------------------------------
cmd_flank <- c(1000, 1000)
cz1_flank <- c(1000, 1000)
cz2_flank <- c(1000, 1000)

# Add C-region highlighting to stacked map
htrack <- 
  HighlightTrack(
    trackList = list(sublocus_track, vh_track, dh_track, jh_track, cm_track, cd_track, cz_track),
    start = c(min(cz1_tab_summary$start, cz1_tab_summary$end)-cmd_flank[1],
              min(cz2_tab_summary$start, cz2_tab_summary$end)-cz1_flank[1],
              min(cmd_tab_summary$start, cmd_tab_summary$end)-cz2_flank[1]),
    end = c(max(cz1_tab_summary$start, cz1_tab_summary$end)+cmd_flank[2],
            max(cz2_tab_summary$start, cz2_tab_summary$end)+cz1_flank[2],
            max(cmd_tab_summary$start, cmd_tab_summary$end)+cz2_flank[2]),
    chromosome = chromosome_id,
    from = 1,
    to = locus_width
    )

# Combine C-region tracks
ctrack_md <-
  HighlightTrack(
    trackList = list(xtrack2, cmd_summary_track, cmd_grouped_track),
    chromosome = "chr6",
    from = min(cmd_tab_summary$start, cmd_tab_summary$end)-cmd_flank[1],
    to = max(cmd_tab_summary$start, cmd_tab_summary$end)+cmd_flank[2]
  )
ctrack_z1 <-
  HighlightTrack(
    trackList = list(xtrack2, cz1_summary_track, cz1_grouped_track),
    chromosome = "chr6",
    from = min(cz1_tab_summary$start, cz1_tab_summary$end)-cz1_flank[1],
    to = max(cz1_tab_summary$start, cz1_tab_summary$end)+cz1_flank[2]
  )
ctrack_z2 <-
  HighlightTrack(
    trackList = list(xtrack2, cz2_summary_track, cz2_grouped_track),
    chromosome = "chr6",
    from = min(cz2_tab_summary$start, cz2_tab_summary$end)-cz2_flank[1],
    to = max(cz2_tab_summary$start, cz2_tab_summary$end)+cz2_flank[2]
  )


# Stack plots together
width_md <- max(cmd_tab_summary$start, cmd_tab_summary$end)+ cmd_flank[2] - 
  (min(cmd_tab_summary$start, cmd_tab_summary$end)- cmd_flank[1])
width_z1 <- max(cz1_tab_summary$start, cz1_tab_summary$end)+ cz1_flank[2] - 
  (min(cz1_tab_summary$start, cz1_tab_summary$end)- cz1_flank[1])
width_z2 <- max(cz2_tab_summary$start, cz2_tab_summary$end)+ cz2_flank[2] - 
  (min(cz2_tab_summary$start, cz2_tab_summary$end)- cz2_flank[1])

width_ratio <- c(width_z2/width_z1, width_md/width_z1)
width_partition <- c(1, width_ratio)
# Stack plots together
plot_unit = "cm"
plot_height <- 24.3*1.5
plot_width <- 20.9*1.5
margin_width <- c(1,1)
fig_height_ratios <- c(1,6,4)
fig_heights <- sapply(fig_height_ratios, 
                      function(x) x/sum(fig_height_ratios) * plot_height)
fig_widths <- sapply(width_partition,
                     function(x) x/sum(width_partition) * plot_width)
map_layout <- grid.layout(
  ncol = 4,
  nrow = 4,
  heights = unit(c(margin_width[2],fig_heights), plot_unit),
  widths = unit(c(margin_width[1],fig_widths), plot_unit)
)

vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
# Add ideogram to first row
pushViewport(viewport(layout.pos.col = 2:4 , layout.pos.row = 2))
plotTracks(itrack, from=locus_start, to=locus_end, add=TRUE)
popViewport(1)
# Add stacked map to second row
pushViewport(viewport(layout.pos.col = 2:4, layout.pos.row = 3))
plotTracks(list(xtrack1, htrack), add=TRUE)
popViewport(1)
# Add C-region maps to third row
pushViewport(viewport(layout.pos.col = 4, layout.pos.row = 4))
plotTracks(ctrack_md, add=TRUE)
popViewport(1)
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 4))
z1_ticks_from <- floor(ctrack_z1@dp@pars$from * 10^-3)/10^-3
z1_ticks_to <- ceiling(ctrack_z1@dp@pars$to * 10^-3)/10^-3
plotTracks(ctrack_z1, add=TRUE, ticksAt = seq(z1_ticks_from,z1_ticks_to,2000))
popViewport(1)
pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 4))
z2_ticks_from <- floor(ctrack_z2@dp@pars$from * 10^-3)/10^-3
z2_ticks_to <- ceiling(ctrack_z2@dp@pars$to * 10^-3)/10^-3
plotTracks(ctrack_z2, add=TRUE, ticksAt = seq(z2_ticks_from,z2_ticks_to,2000))
popViewport(1)
# Add subfigure labels
gp_label <- gpar(fontfamily = titlefont, fontsize = fontsize_base * 3)
pushViewport(viewport(layout.pos.col = 1:length(map_layout$widths),
                      layout.pos.row = 1:length(map_layout$heights)))
grid.text(label = c("A", "B", "C"), x = 0.6, 
          y = c(plot_height + 1, 
                plot_height - sum(fig_heights[1]) + 0,
                plot_height - sum(fig_heights[1:2])) - 1,
          just = "centre", default.units = "cm",
          gp = gp_label)
popViewport(1)

#------------------------------------------------------------------------------
# ADD ZOOM ANNOTATIONS
#------------------------------------------------------------------------------

# Define graphical and plotting settings
gpoly <- function(x1, x2, y1, y2, col, alpha=0.1,
                  lwd = 1){
  gp_fill <- gpar(col=col, fill=col, alpha=alpha,
                  lty=2, lwd=lwd, linejoin="round")
  gp_unfill <- gpar(col=col, lty=2, lwd=lwd, linejoin="round")
  grid.polygon(x = c(x1, rev(x2), x1[1]),
               y = c(y1, y2, y1[1]),
               default.units = "cm", gp = gp_unfill)
  grid.polygon(x = c(x1, rev(x2), x1[1]),
               y = c(y1, y2, y1[1]),
               default.units = "cm", gp = gp_fill)
}

# Define zoom regions
lower_track_y <- rep(12.20, 2)
lower_hilight_y <- rep(13.42, 2)
cz1_track_x <- c(2.230, 8.22)
cz1_hilight_x <- c(2.162, 2.961)
cz2_track_x <- c(10.05, 14.88)
cz2_hilight_x <- c(27.79, 28.45)
cmd_track_x <- c(16.585, 31.935)
cmd_hilight_x <- c(30.08, 31.84)

higher_track_y <- lower_track_y + 19.8
higher_hilight_y <- lower_hilight_y + 20.15
higher_mid_y <- (higher_track_y + 2*higher_hilight_y)/3
higher_track_x <- c(2.24, 31.925)
higher_hilight_x <- c(21.648, 21.949)
higher_mid_x <- rep(mean(higher_hilight_x), 2)

# Add lower zoom lines
pushViewport(viewport(layout.pos.col = 1:length(map_layout$widths),
                      layout.pos.row = 1:length(map_layout$heights)))
gpoly(cz1_track_x, cz1_hilight_x, lower_track_y, lower_hilight_y,
      colours[["HL2"]])
gpoly(cz2_track_x, cz2_hilight_x, lower_track_y, lower_hilight_y,
      colours[["HL2"]])
gpoly(cmd_track_x, cmd_hilight_x, lower_track_y, lower_hilight_y,
      colours[["HL2"]])

# Add higher zoom lines, with reversal indicator
gpoly(higher_track_x, higher_mid_x, higher_track_y, higher_mid_y,
      colours[["HL1"]])
gpoly(higher_mid_x, higher_hilight_x, higher_mid_y, higher_hilight_y,
      colours[["HL1"]])
grid.lines(x = higher_hilight_x + c(-0.1, 0.1),
          y = higher_hilight_y + 2.65,
          arrow = arrow(ends="first", length=unit(0.2, "cm"),
                        type="closed"),
          gp = gpar(col=colours[["HL1"]], fill=colours[["HL1"]], lwd = 2),
          default.units = "cm")
popViewport(1)

#------------------------------------------------------------------------------
# ADD TEXT ANNOTATIONS
#------------------------------------------------------------------------------

# Define graphics settings
label_names <- c("VH", "DH", "JH", "CZ", "CM", "CD")
gp_txt_list <- lapply(label_names, function(n)
  gpar(fontfamily = titlefont, fontsize = fontsize_base * 1.7,
       col = colours[[n]], lwd=1)
)
names(gp_txt_list) <- label_names

# Define label expressions
cz1_labels <- c(paste0("'C'[zeta]*'", seq(4), "'"), "TM1", "TM2")
cz2_labels <- c(paste0("'C'[zeta]*'", seq(4), "'"), "TM1", "TM2")
cmd_labels <- c(paste0("'C'[mu]*'", seq(4), "'"), "TM1", "TM2",
                "'C'[delta]*'1'", paste0("'C'[delta]*'", seq(2,4), "a'"),
                paste0("'C'[delta]*'", seq(2,4), "b'"),
                paste0("'C'[delta]*'", seq(5,7), "'"), "TM1", "TM2")

# Define x and y co-ordinates
halfrep <- function(e1, e2, n){
  v <- rep(c(e1, e2), n %/% 2)
  if (n %% 2 == 0) return(v)
  if (n %% 2 == 1) return(c(v, e1))
  stop("Not an integer.")
}
label_space <- 1.3
label_start <- 1.8
label_upper_y <- 10.5
label_lower_y <- 6.6
label_cz1_x <- rev(8.4 - rep(label_start + 0:length(cz1_labels) * label_space, 
                   each=2)[1:length(cz1_labels)])
label_cz1_y <- halfrep(label_upper_y, label_lower_y, length(cz1_labels))
label_cz2_x <- rev(7.5 - rep(label_start + 0:length(cz2_labels) * label_space, 
                   each=2)[1:length(cz2_labels)])
label_cz2_y <- halfrep(label_upper_y, label_lower_y, length(cz2_labels))
label_cmd_x <- rev(17.7 - rep(label_start + 0:length(cmd_labels) * label_space*1.2, 
                   each=2)[1:length(cmd_labels)])
label_cmd_y <- halfrep(label_upper_y, label_lower_y, length(cmd_labels))

gtext <- function(lab, x, y, tm_gp, just=c(0.5, 1)){
  if (lab %in% c("VH", "DH", "JH")){ gp <- gp_txt_list[[lab]]
  } else if (grepl("zeta", lab)){ gp <- gp_txt_list[["CZ"]]
  } else if (grepl("mu", lab)){ gp <- gp_txt_list[["CM"]]
  } else if (grepl("delta", lab)){ gp <- gp_txt_list[["CD"]]
  } else if (grepl("TM", lab)){ gp <- tm_gp
  } else {stop("Unrecognised label type.")}
  grid.text(label = bquote(.(parse(text=lab))), x = x, y = y, gp = gp,
            just = just, default.units = "cm")
}

# Print labels
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 4))
for (n in 1:length(cz1_labels)){
  lab <- cz1_labels[n]
  gtext(lab, label_cz1_x[n], label_cz1_y[n], gp_txt_list[["CZ"]])
}
popViewport(1)
pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 4))
for (n in 1:length(cz2_labels)){
  lab <- cz2_labels[n]
  gtext(lab, label_cz2_x[n], label_cz2_y[n], gp_txt_list[["CZ"]])
}
popViewport(1)
pushViewport(viewport(layout.pos.col = 4, layout.pos.row = 4))
tm_gp <- gp_txt_list[["CM"]]
for (n in 1:length(cmd_labels)){
  lab <- cmd_labels[n]
  if (grepl("delta", lab)) tm_gp <- gp_txt_list[["CD"]]
  gtext(lab, label_cmd_x[n], label_cmd_y[n], tm_gp)
}
popViewport(1)

# Add line markers

gline <- function(lab, x1, x2, y1, y2, tm_gp){
  if (lab %in% c("VH", "DH", "JH")){ gp <- gp_txt_list[[lab]]
  } else if (grepl("zeta", lab)){ gp <- gp_txt_list[["CZ"]]
  } else if (grepl("mu", lab)){ gp <- gp_txt_list[["CM"]]
  } else if (grepl("delta", lab)){ gp <- gp_txt_list[["CD"]]
  } else if (grepl("TM", lab)){ gp <- tm_gp
  } else {stop("Unrecognised label type.")}
  grid.lines(x = c(x1, x2), y = c(y1, y2), gp = gp, default.units = "cm")
}
y1_down <- -0.7
y1_up <- 0.1
y2_down <- -0.98
y2_up <- 0.55

pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 4))
y1 <- label_cz1_y + halfrep(y1_down, y1_up, length(cz1_labels))
y2 <- label_cz1_y + halfrep(y2_down, y2_up, length(cz1_labels))
x1 <- label_cz1_x
x2 <- x1 + c(.59, .92, -.06, .48, .23, .54)
for (n in 1:length(cz1_labels)){
  lab <- cz1_labels[n]
  gline(lab, x1[n], x2[n], y1[n], y2[n], gp_txt_list[["CZ"]])
}
popViewport(1)

pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 4))
y1 <- label_cz2_y + halfrep(y1_down, y1_up, length(cz2_labels))
y2 <- label_cz2_y + halfrep(y2_down, y2_up, length(cz2_labels))
x1 <- label_cz2_x
x2 <- x1 + c(1.32, 1.60, 0.70, 1.08, .27, .56)
for (n in 1:length(cz2_labels)){
  lab <- cz2_labels[n]
  gline(lab, x1[n], x2[n], y1[n], y2[n], gp_txt_list[["CZ"]])
}
popViewport(1)

pushViewport(viewport(layout.pos.col = 4, layout.pos.row = 4))
y1 <- label_cmd_y + halfrep(y1_down, y1_up, length(cmd_labels))
y2 <- label_cmd_y + halfrep(y2_down, y2_up, length(cmd_labels))
x1 <- label_cmd_x
x2 <- x1 + c(2.86,3.92,2.78,3.28,2.54,4.71,
             3.66,4.02,2.92,3.24,
             3.12, 3.50, 2.29, 2.61, 1.45, 1.77, .53, .69)
tm_gp <- gp_txt_list[["CM"]]
for (n in 1:length(cmd_labels)){
  lab <- cmd_labels[n]
  if (grepl("delta", lab)) tm_gp <- gp_txt_list[["CD"]]
  gline(lab, x1[n], x2[n], y1[n], y2[n], tm_gp)
}
popViewport(1)

#------------------------------------------------------------------------------
# SAVE FIGURE
#------------------------------------------------------------------------------

plt <- grid.grab()

savefig(plot = plt, filename = filename, 
        height = plot_height + margin_width[2], 
        width = plot_width + margin_width[1])

#------------------------------------------------------------------------------
# FIX GVIZ INSTALLATION
#------------------------------------------------------------------------------
BiocManager::install("Gviz")
if (packageVersion("Gviz") != gvizV1) stop("Error: failed to restore Gviz installation.")