###############################################################################
## AUX FILE                                                                  ##
## Configure Gviz display defaults                                           ##
###############################################################################

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")

# Specify basic options
options(ucscChromosomeNames=FALSE)

# Create new display scheme
scheme <- getScheme()

# Gviz-specific font settings
fontcolor_ideogram <- "#606060"
fontsize_base <- 10
fontscale_ideogram <- 1.2
fontscale_gaxis <- 1.0
fontscale_annotation_title <- 1.2
fontscale_annotation_group <- 1.0

# IdeogramTrack
scheme$IdeogramTrack$fontfamily = font
scheme$IdeogramTrack$cex = fontscale_ideogram
scheme$IdeogramTrack$bevel = 1
scheme$IdeogramTrack$fontcolor = fontcolour_mid
scheme$IdeogramTrack$fontsize = fontsize_base
scheme$IdeogramTrack$col = paste0(colours[["HL1"]], "FF")
scheme$IdeogramTrack$fill = paste0(colours[["HL1"]], "20")
scheme$IdeogramTrack$showId = TRUE

# GenomeAxisTrack
scheme$GenomeAxisTrack$add53 <- TRUE
scheme$GenomeAxisTrack$add35 <- TRUE
scheme$GenomeAxisTrack$labelPos <- "below"
scheme$GenomeAxisTrack$lwd <- 3
scheme$GenomeAxisTrack$cex <- fontscale_gaxis
scheme$GenomeAxisTrack$fontsize <- fontsize_base
scheme$GenomeAxisTrack$fontfamily <- font
scheme$GenomeAxisTrack$fontfamily.title <- titlefont
scheme$GenomeAxisTrack$distFromAxis = 3
scheme$GenomeAxisTrack$cex = 1
scheme$GenomeAxisTrack$size=0.9

# AnnotationTrack
scheme$AnnotationTrack$col <- NULL
scheme$AnnotationTrack$fontfamily.title <- titlefont
scheme$AnnotationTrack$fontfamily <- font
scheme$AnnotationTrack$cex.title <- fontscale_annotation_title
scheme$AnnotationTrack$stacking <- "dense"
scheme$AnnotationTrack$VH <- colours[["VH"]]
scheme$AnnotationTrack$DH <- colours[["DH"]]
scheme$AnnotationTrack$JH <- colours[["JH"]]
scheme$AnnotationTrack$CD <- colours[["CD"]]
scheme$AnnotationTrack$CM <- colours[["CM"]]
scheme$AnnotationTrack$CZ <- colours[["CZ"]]
scheme$AnnotationTrack$Sublocus <- colours[["GR"]]
scheme$AnnotationTrack$background.title = colours[["GR"]]
scheme$AnnotationTrack$background.panel = paste0(colours[["GR"]], "20")
scheme$AnnotationTrack$shape <- "box"
scheme$AnnotationTrack$cex.group = 1

# GeneRegionTrack
scheme$GeneRegionTrack$stacking <- "squish"
scheme$GenomeAxisTrack$fontfamily <- font
scheme$GenomeAxisTrack$fontfamily.title <- titlefont
scheme$GeneRegionTrack$background.title = paste0(colours[["GR"]], "80")
scheme$GeneRegionTrack$background.panel = paste0(colours[["GR"]], "10")
scheme$GeneRegionTrack$col <- NULL

# HighlightTrack
scheme$GenomeAxisTrack$fontfamily <- font
scheme$GenomeAxisTrack$fontfamily.title <- titlefont
scheme$HighlightTrack$col <- paste0(colours[["HL2"]], "FF")
scheme$HighlightTrack$fill <- paste0(colours[["HL2"]], "20")
scheme$HighlightTrack$inBackground <- TRUE

# AlignmentsTrack
scheme$GenomeAxisTrack$fontfamily <- font
scheme$GenomeAxisTrack$fontfamily.title <- titlefont
scheme$AlignmentsTrack$background.title = colours[["GR"]]
scheme$AlignmentsTrack$background.panel = colours[["WH"]]


# Activate scheme
addScheme(scheme, "locus_map_scheme")
options(Gviz.scheme = "locus_map_scheme")
