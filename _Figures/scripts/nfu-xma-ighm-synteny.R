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
locus_path_nfu <- "../_Data/locus/complete/nfu.fasta"
locus_path_xma <- "../_Data/locus/complete/xma.fasta"

ch_ranges_path_nfu <- "../_Data/ranges/nfu/nfu_ch_ranges.tsv"
ch_ranges_path_xma <- "../_Data/ranges/xma/xma_ch_ranges.tsv"

# Configure output
filename_base <- "nfu-xma-ighm-synteny"

#------------------------------------------------------------------------------
# PREPARE CM-REGION RANGES
#------------------------------------------------------------------------------

if (!exists("igh_nfu")){
  igh_nfu <- readDNAStringSet(locus_path_nfu)
}
if (!exists("igh_xma")){
  igh_xma <- readDNAStringSet(locus_path_xma)
}
igh_loci <- c(igh_nfu, igh_xma)

# Set sublocus information
subloci_labels <- c("IGH1", "IGH2")
subloci_seqnames <- c("locus", "locus")
subloci_strands <- c("+", "-")

# Get CM ranges for annotation (excluding TM2)
cm_tab_nfu <- suppressMessages(read_tsv(ch_ranges_path_nfu)) %>%
  filter(grepl("IGH1M", label), !grepl("TM2", label))
cm_tab_xma <- suppressMessages(read_tsv(ch_ranges_path_xma)) %>%
  filter(grepl("IGHM", label), !grepl("TM2", label))
cm_tab <- bind_rows(cm_tab_nfu %>% mutate(species = "nfu"),
                    cm_tab_xma %>% mutate(species = "xma"))

# Infer CM-region boundaries
flank_width <- 500
cm_regions_tab <- cm_tab %>% group_by(species, seqnames, strand) %>%
  summarise(start = min(start)-flank_width, end = max(end)+flank_width,
            width = end - start + 1) %>%
  mutate(labels = paste0("igh_", species))

# Convert to ranges and extract CM-regions from locus
cm_regions_ranges <- cm_regions_tab %>% as("GRanges")
cm_regions <- getSeq(igh_loci, cm_regions_ranges)
names(cm_regions) <- cm_regions_tab$labels

#------------------------------------------------------------------------------
# CREATE DECIPHER DB AND CALCULATE SYNTENY
#------------------------------------------------------------------------------

DBPath <- tempfile()

DBConn <- dbConnect(SQLite(),
                    DBPath)

Seqs2DB(seqs = cm_regions[1], type = "XStringSet", dbFile = DBConn,
        identifier = names(cm_regions)[1], verbose = FALSE)
Seqs2DB(seqs = cm_regions[2], type = "XStringSet", dbFile = DBConn,
        identifier = names(cm_regions)[2], verbose = FALSE)
dbDisconnect(DBConn)

SyntenyObject <- FindSynteny(dbFile = DBPath,
                             verbose = FALSE)

#------------------------------------------------------------------------------
# EXTRACT SYNTENY INFORMATION AND MAKE ALIGNMENT PLOT
#------------------------------------------------------------------------------

SyntenyTableA <- SyntenyObject[2][[1]] %>% as.data.frame %>% as.tbl
SyntenyTableB <- SyntenyObject[3][[1]] %>% as.data.frame %>% as.tbl %>%
  mutate(end1 = ifelse(strand == 0, start1 + width - 1, start1 - width + 1),
         end2 = ifelse(strand == 0, start2 + width - 1, start2 - width + 1))
nfu_start <- cm_regions_tab %>% filter(species == "nfu") %>% pull(start)
xma_start <- cm_regions_tab %>% filter(species == "xma") %>% pull(start)
nfu_end <- cm_regions_tab %>% filter(species == "nfu") %>% pull(end)
xma_end <- cm_regions_tab %>% filter(species == "xma") %>% pull(end)

# Create synteny plot with region labels
highlight_offset <- 1000
g <- ggplot() + 
  geom_segment(data = SyntenyTableB %>% filter(strand == 0),
               aes(x = start1, xend = end1, y = start2, yend = end2),
               colour = "black") +
  geom_segment(data = SyntenyTableB %>% filter(strand == 1),
               aes(x = start1, xend = end1, y = end2, yend = start2),
               colour = "red") +
  geom_rect(data = cm_tab_nfu,
            aes(ymin = 0, ymax = Inf, 
                xmin = start - nfu_start, xmax = end - nfu_start, 
                fill = colours[["CM"]]), 
                alpha = 0.3) +
  geom_rect(data = cm_tab_xma,
            aes(xmin = 0, xmax = Inf, 
                ymin = start - xma_start, ymax = end - xma_start, 
                fill = colours[["CM"]]), 
            alpha = 0.3)
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

plot_height <- 20
plot_ratio <- 1.4

savefig(g, filename, height = plot_height, ratio = plot_ratio)