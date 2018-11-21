###############################################################################
## FIGURE                                                                    ##
## Nothobranchius furzeri IgH1/2 dotplot                                     ##
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
locus_path <- "../_Data/locus/complete/nfu.fasta"
vh_ranges_path <- "../_Data/ranges/nfu/nfu_vh_ranges.tsv"
dh_ranges_path <- "../_Data/ranges/nfu/nfu_dh_ranges.tsv"
jh_ranges_path <- "../_Data/ranges/nfu/nfu_jh_ranges.tsv"
ch_ranges_path <- "../_Data/ranges/nfu/nfu_ch_ranges.tsv"

# Configure output
filename <- "nfu-locus-dots"

#------------------------------------------------------------------------------
# IDENTIFY SUBLOCUS RANGES
#------------------------------------------------------------------------------

if (!exists("igh_nfu")){
  igh_nfu <- readDNAStringSet(locus_path)
}

# Set sublocus information
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

# Extend subloci into flanking regions
flank_width <- 1000
subloci <- subloci %>% 
  mutate(start = pmax(0, start - flank_width),
         end = pmin(end + flank_width, width(igh_nfu)),
         width = end - start + 1)

# Convert to ranges and extract sequences from locus
sublocus_ranges <- subloci %>% as("GRanges")
sublocus_seqs <- getSeq(igh_nfu, sublocus_ranges)

#------------------------------------------------------------------------------
# CREATE DECIPHER DB AND CALCULATE SYNTENY
#------------------------------------------------------------------------------

DBPath <- tempfile()

DBConn <- dbConnect(SQLite(),
                    DBPath)

Seqs2DB(seqs = sublocus_seqs[1], type = "XStringSet", dbFile = DBConn,
        identifier = "IGH1", verbose = FALSE)
Seqs2DB(seqs = sublocus_seqs[2], type = "XStringSet", dbFile = DBConn,
        identifier = "IGH2", verbose = FALSE)
dbDisconnect(DBConn)

SyntenyObject <- FindSynteny(dbFile = DBPath,
                             verbose = FALSE)

#------------------------------------------------------------------------------
# EXTRACT SYNTENY INFORMATION AND DETERMINE HIGHLIGHTING REGIONS
#------------------------------------------------------------------------------

SyntenyTable <- SyntenyObject[2][[1]] %>% as.data.frame %>% as.tbl
igh1_start <- subloci %>% filter(label == "IGH1") %>% pull(start)
igh2_start <- subloci %>% filter(label == "IGH2") %>% pull(start)
igh1_end <- subloci %>% filter(label == "IGH1") %>% pull(end)
igh2_end <- subloci %>% filter(label == "IGH2") %>% pull(end)

# Create V/D/J/C-region tables
make_region_table <- function(tab, region, exclude, class){
  out_tab <- tab %>%
    filter(grepl(region, label), !(label %in% exclude)) %>%
    summarise(label = paste(region, class, sep = "-"),
              start = min(start), end = max(end),
              strand = first(strand), seqnames = first(seqnames)) %>%
    mutate(width = end - start + 1)
  exc_tab <- tab %>%
    filter(grepl(region, label), label %in% exclude)
  return(bind_rows(out_tab, exc_tab) %>% mutate(class = class))
}

igh1_regions <- bind_rows(
  make_region_table(vh_tab, "IGH1", "IGH1V1-07", "VH"),
  make_region_table(dh_tab, "IGH1", "IGH1D01", "DH"),
  make_region_table(jh_tab, "IGH1", "IGH1J01", "JH"),
  make_region_table(cm_tab, "IGH1", "", "CM"),
  make_region_table(cd_tab, "IGH1", "", "CD")
) %>%
  mutate(class = factor(class, levels = names(colours)),
         start = start - igh1_start, end = end - igh1_start)
igh2_regions <- bind_rows(
  make_region_table(vh_tab, "IGH2", "", "VH"),
  make_region_table(dh_tab, "IGH2", "", "DH"),
  make_region_table(jh_tab, "IGH2", "", "JH"),
  make_region_table(cm_tab, "IGH2", "", "CM"),
  make_region_table(cd_tab, "IGH2", "", "CD")
) %>%
  mutate(class = factor(class, levels = names(colours)),
         start = igh2_end - start, end = igh2_end - end)

# Determine highlighting areas for main regions
igh1_r <- igh1_regions %>% filter(grepl("IGH1-", label)) %>% 
  select(-label, -strand, -width)
igh2_r <- igh2_regions %>% filter(grepl("IGH2-", label)) %>% 
  select(-label, -strand, -width)
igh_highlight <- full_join(igh1_r, igh2_r, by=c("class", "seqnames")) %>%
  select(seqnames, class, everything())


# Create synteny plot with region labels
g <- ggplot() + 
  geom_segment(data = SyntenyTable %>% filter(strand == 0),
               aes(x = start1, xend = end1, y = start2, yend = end2),
               colour = "black") +
  geom_segment(data = SyntenyTable %>% filter(strand == 1),
               aes(x = start1, xend = end1, y = end2, yend = start2),
               colour = "red") +
  geom_rect(data = igh1_regions,
            aes(ymin = -4000, ymax = -2000, xmin = start, xmax = end,
                colour = class, fill = class), alpha = 0.5) +
  geom_rect(data = igh2_regions,
            aes(xmin = -4000, xmax = -2000, ymin = start, ymax = end,
                colour = class, fill = class), alpha = 0.5) +
  geom_rect(data = igh_highlight,
            aes(xmin = -2000, xmax = end.x, ymin = start.y, ymax = end.y,
                fill = class), alpha = 0.3) +
  geom_rect(data = igh_highlight,
            aes(xmin = start.x, xmax = end.x, ymin = -2000, ymax = end.y,
                fill = class), alpha = 0.3) +
  geom_rect(data = igh_highlight,
            aes(xmin = start.x, xmax = end.x, ymin = start.y, ymax = end.y,
                fill = class), alpha = 0.3) +
  scale_colour_manual(values = as.character(colours), name = "Region") + 
  scale_fill_manual(values = as.character(colours), name = "Region") +
  coord_fixed() + theme_minimal() +
  xlab("IGH1 Co-ordinate") + ylab("IGH2 Co-ordinate") +
  theme(
    legend.position = "top",
    axis.text = element_text(size = fontsize_base, family = "CM Sans", 
                             colour = "black"),
    axis.title.y = element_text(size = fontsize_base * fontscale_title,
                              family = font, colour = "black",
                              margin = margin(r=2,unit="mm")),
    axis.title.x = element_text(size = fontsize_base * fontscale_title,
                                family = font, colour = "black",
                                margin = margin(t=5,unit="mm")),
    legend.text = element_text(size = fontsize_base, 
                               family = font, colour = "black"),
    legend.title = element_text(size = fontsize_base,
                                family = font, colour = "black",
                                face = "bold")
  )

# TODO: Add highlighted regions (currently done manually)

#------------------------------------------------------------------------------
# SAVE FIGURE
#------------------------------------------------------------------------------

plot_height <- 15
plot_ratio <- 1.3

savefig(g, filename, height = plot_height, ratio = plot_ratio)