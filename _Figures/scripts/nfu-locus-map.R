###############################################################################
## FIGURE                                                                    ##
## Nothobranchius furzeri IgH locus map                                      ##
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
genome_path <- "/home/will/Documents/data/genome/NFZv2.0.fasta"
locus_path <- "/home/will/Documents/data/igh/locus/locus.fa"
vh_ranges_path <- "/home/will/Documents/data/igh/locus/vsearch/nfu_vh_ranges.tsv"
dh_ranges_path <- "/home/will/Documents/data/igh/locus/dsearch/nfu_dh_ranges.tsv"
jh_ranges_path <- "/home/will/Documents/data/igh/locus/jsearch/nfu_jh_ranges.tsv"
ch_ranges_path <- "/home/will/Documents/data/igh/locus/csearch/nfu_ch_ranges.tsv"

# Configure output
filename <- "nfu-locus-map"

#------------------------------------------------------------------------------
# LOCUS IDEOGRAM
#------------------------------------------------------------------------------

# Specify locus and chromosome information
genome_name <- "NFZ2.0"
chromosome_id <- "chr6"
chromosome_name <- "chr6"
chromosome_length <- 61955777
locus_start <- 43205521 + 98762 - 1 # (on chromosome)
locus_end <- 43205521 + 253166 - 1

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
locus_width <- 306408
subloci_labels <- c("IGH1", "IGH2")
subloci_seqnames <- c("locus", "locus")
subloci_strands <- c("+", "-")

# Get V/D/J/C ranges for annotation
vh_tab <- suppressMessages(read_tsv(vh_ranges_path))
ch_tab <- suppressMessages(read_tsv(ch_ranges_path))
jh_tab <- suppressMessages(read_tsv(jh_ranges_path))
dh_tab <- suppressMessages(read_tsv(dh_ranges_path))

# Split up CH ranges
cd_tab <- ch_tab %>% filter(grepl("IGH.D", label))
cm_tab <- ch_tab %>% filter(grepl("IGH.M", label))
cz_tab <- ch_tab %>% filter(grepl("IGH.Z", label))

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
# ZOOMED C-REGION MAP
#------------------------------------------------------------------------------

# Re-annotate C-range table and split into subloci
ch_tab <- ch_tab %>%
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

# Add JH ranges
# Re-annotate range table and split into subloci
jh_tab <- jh_tab %>% mutate(feature = "JH",
                            exon = paste0("JH-", sub("IGH\\dJ", "", label)))
cj_tab <- bind_rows(ch_tab, jh_tab)
cj_tab_subloci <- lapply(subloci_labels, function(l) 
  cj_tab %>% filter(grepl(l, label)) %>% 
    filter(feature != "JH" | pmin(abs(min(ch_tab_subloci[[l]]$start) - start),
                                  abs(start - max(ch_tab_subloci[[l]]$start))) < 10000))
names(cj_tab_subloci) <- subloci_labels

# Prepare tracks
# Make overall and transcript-specific tracks
sublocus_summary_track <- lapply(subloci_labels, function(l) 
  AnnotationTrack(cj_tab_subloci[[l]],
                  feature = cj_tab_subloci[[l]]$feature,
                  group = cj_tab_subloci[[l]]$exon,
                  name = "Exons",
                  groupAnnotation = "group",
                  chromosome = chromosome_id,
                  stackHeight=0.4,
                  size = 3) # Give enough height to manually add labels later
)
names(sublocus_summary_track) <- subloci_labels

sublocus_grouped_track <- lapply(subloci_labels, function(l)
  GeneRegionTrack(ch_tab_subloci_grouped[[l]],
                  group = ch_tab_subloci_grouped[[l]]$transcript,
                  feature = ch_tab_subloci_grouped[[l]]$feature,
                  transcript = ch_tab_subloci_grouped[[l]]$transcript,
                  name = "Isoforms",
                  transcriptAnnotation = "transcript",
                  chromosome = chromosome_id,
                  stacking = "squish")
)
names(sublocus_grouped_track) <- subloci_labels

#------------------------------------------------------------------------------
# COMBINE AND PLOT TRACKS
#------------------------------------------------------------------------------

# Add C-region highlighting to stacked map
htrack <- lapply(subloci_labels, function(l)
  HighlightTrack(
    trackList = list(sublocus_track, vh_track, dh_track, jh_track, cm_track, cd_track, cz_track),
    start = min(cj_tab_subloci[[l]]$start, cj_tab_subloci[[l]]$end)-2000,
    end = max(cj_tab_subloci[[l]]$start, cj_tab_subloci[[l]]$end)+2000,
    chromosome = "chr6",
    from = 1,
    to = locus_width
    )
)
names(htrack) <- subloci_labels

# Combine C-region tracks
ctrack <- lapply(subloci_labels, function(l)
  HighlightTrack(
    trackList = list(xtrack2, sublocus_summary_track[[l]], sublocus_grouped_track[[l]]),
    chromosome = "chr6",
    from = min(cj_tab_subloci[[l]]$start, cj_tab_subloci[[l]]$end)-2000,
    to = max(cj_tab_subloci[[l]]$start, cj_tab_subloci[[l]]$end)+2000
  ))
names(ctrack) <- subloci_labels

# Stack plots together
plot_unit = "cm"
plot_height <- 24.3*1.5
plot_width <- 20.9*1.5
margin_width <- c(1,1)
fig_height_ratios <- c(1,6,4)
fig_heights <- sapply(fig_height_ratios, 
                           function(x) x/sum(fig_height_ratios) * plot_height)
map_layout <- grid.layout(
  ncol = 2,
  nrow = 4,
  heights = unit(c(margin_width[2],fig_heights), plot_unit),
  widths = unit(c(margin_width[1],plot_width), plot_unit)
)
vtop <- viewport(layout = map_layout)
grid.newpage()
pushViewport(vtop)
# Add ideogram to first row
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
plotTracks(itrack, from=locus_start, to=locus_end, add=TRUE)
popViewport(1)
# Add stacked map to second row
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3))
plotTracks(list(xtrack1, htrack[["IGH1"]]), add=TRUE)
popViewport(1)
# Add C-region map to third row
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 4))
plotTracks(ctrack[["IGH1"]], add=TRUE)
popViewport(1)
# Add subfigure labels
gp <- gpar(fontfamily = titlefont, fontsize = fontsize_base * 3)
pushViewport(viewport(layout.pos.col = 1:length(map_layout$widths),
                      layout.pos.row = 1:length(map_layout$heights)))
grid.text(label = c("A", "B", "C"), x = 0.6, 
          y = c(plot_height + 1, 
                plot_height - sum(fig_heights[1]) + 0,
                plot_height - sum(fig_heights[1:2])) - 1,
          just = "centre", default.units = "cm",
          gp = gp)
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
cmd_track_x <- c(2.26, 31.935)
cmd_hilight_x <- c(14.38, 16.74)

higher_track_y <- lower_track_y + 19.8
higher_hilight_y <- lower_hilight_y + 20.15
higher_track_x <- c(2.24, 31.925)
higher_hilight_x <- c(23.88, 24.02)

# Add lower zoom lines
pushViewport(viewport(layout.pos.col = 1:length(map_layout$widths),
                      layout.pos.row = 1:length(map_layout$heights)))
gpoly(cmd_track_x, cmd_hilight_x, lower_track_y, lower_hilight_y,
      colours[["HL2"]])

# Add higher zoom lines
gpoly(higher_track_x, higher_hilight_x, higher_track_y, higher_hilight_y,
      colours[["HL1"]])
grid.lines(x = higher_hilight_x + c(-0.1, 0.1),
           y = higher_hilight_y + 2.65,
           arrow = arrow(ends="last", length=unit(0.2, "cm"),
                         type="closed"),
           gp = gpar(col=colours[["HL1"]], fill=colours[["HL1"]], lwd = 2),
           default.units = "cm")
popViewport(1)

#------------------------------------------------------------------------------
# ADD TEXT ANNOTATIONS
#------------------------------------------------------------------------------

# Define graphics settings
label_names <- c("VH", "DH", "JH", "CM", "CD")
gp_txt_list <- lapply(label_names, function(n)
  gpar(fontfamily = titlefont, fontsize = fontsize_base * 1.7,
       col = colours[[n]], lwd=1)
)
names(gp_txt_list) <- label_names

# Define label expressions
cm_labels <- c(paste0("'C'[mu]*'", seq(4), "'"), "TM1")
cd_labels <- c("TM2", "'C'[delta]*'1'", 
               paste0("'C'[delta]*'", seq(2,4), "a'"),
                paste0("'C'[delta]*'", seq(2,4), "b'"),
                paste0("'C'[delta]*'", seq(5,7), "'"), "TM1", "TM2")

# Define x and y co-ordinates
halfrep <- function(e1, e2, n){
  v <- rep(c(e1, e2), n %/% 2)
  if (n %% 2 == 0) return(v)
  if (n %% 2 == 1) return(c(v, e1))
  stop("Not an integer.")
}

label_upper_y <- 10.5
label_lower_y <- 6.6
label_cm_x <- rev(14.5 - rep(2.8 + 0:length(cm_labels) * 2.3, each=2)[1:length(cm_labels)])
label_cd_x <- rev(34 - rep(2.8 + 0:length(cd_labels) * 2.3, each=2)[1:length(cd_labels)])
label_cm_y <- halfrep(label_upper_y, label_lower_y, length(cm_labels))
label_cd_y <- halfrep(label_upper_y, label_lower_y, length(cd_labels))

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
popViewport(1)
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 4))
for (n in 1:length(cm_labels)){
  lab <- cm_labels[n]
  gtext(lab, label_cm_x[n], label_cm_y[n], gp_txt_list[["CM"]])
}
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 4))
tm_gp = gp_txt_list[["CM"]]
for (n in 1:length(cd_labels)){
  lab <- cd_labels[n]
  if (grepl("delta", lab)) tm_gp <- gp_txt_list[["CD"]]
  gtext(lab, label_cd_x[n], label_cd_y[n], tm_gp)
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
y1 <- label_cm_y + halfrep(y1_down, y1_up, length(cm_labels))
y2 <- label_cm_y + halfrep(y2_down, y2_up, length(cm_labels))
x1 <- label_cm_x
x2 <- x1 + c(0.62,-0.11,0.38,-1.32,0.54)
pushViewport(viewport(layout.pos.col = 4, layout.pos.row = 4))
for (n in 1:length(cm_labels)){
  lab <- cm_labels[n]
  gline(lab, x1[n], x2[n], y1[n], y2[n], gp_txt_list[["CM"]])
}

popViewport(1)
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 4))
y1 <- label_cd_y + halfrep(y1_down, y1_up, length(cd_labels))
y2 <- label_cd_y + halfrep(y2_down, y2_up, length(cd_labels))
x1 <- label_cd_x
x2 <- x1 + c(0.85, -0.29, 0.31, -1.49, -1.02,
             3.72, 4.30, 2.49, 
             3.01, 1.35, 1.83, .03, .44)
tm_gp = gp_txt_list[["CM"]]
pushViewport(viewport(layout.pos.col = 4, layout.pos.row = 4))
for (n in 1:length(cd_labels)){
  lab <- cd_labels[n]
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
          