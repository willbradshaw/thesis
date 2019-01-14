###############################################################################
## FIGURE                                                                    ##
## Multispecies C-region maps                                                ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

rm(list = ls())

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")
source("aux/gviz.R")

# Configure species names
species_names_long <- c("Oryzias latipes", "Cyprinodon variegatus",
                        "Fundulus heteroclitus", "Poecilia formosa",
                        "Poecilia reticulata", "Kryptolebias marmoratus",
                        "Austrofundulus limnaeus", "Pachypanchax playfairii",
                        "Callopanchax toddi", "Aphyosemion australe",
                        "Nothobranchius orthonotus")
species_names_short <- sapply(str_split(species_names_long, " "), function(s)
  paste0(substr(s[1], 1, 1), ". ", s[2]))
species_ids <- sapply(str_split(species_names_long, " "), function(s)
  paste0(tolower(substr(s[1], 1, 1)), tolower(substr(s[2], 1, 2))))

# Configure input paths
range_paths <- lapply(species_ids, function(s)
  paste0("../_Data/ranges/", s, "/", s, "-ch-ranges.tsv"))
names(range_paths) <- species_ids

# Configure output
#filename <- "nfu-locus-map"


#------------------------------------------------------------------------------
# STACKED LOCUS MAP
#------------------------------------------------------------------------------

# Get C ranges for annotation
ch_tabs <- lapply(range_paths, function(p) 
  suppressMessages(read_tsv(p)) %>% # Read ranges
    mutate(isotype = sub(".*IGH([MDZ]).*", "C\\1", label), # Identify isotype
           region = sub(".*(IGH.*)-.*", "\\1", label), # Identify constant region
           exon = sub(".*IGH.*-(.*)_?.*", "\\1", label)) %>% # Identify exon
    group_by(seqnames, region, strand, isotype) %>% # Group by scaffold and region
    arrange(start) %>%
    summarise(exon = paste0(exon, collapse=""), start = min(start), 
              end = max(end), width = end - start + 1) %>%
    ungroup() %>% arrange(seqnames, start) %>%
    mutate(exon = gsub("TM1", 8, exon),
           exon = gsub("TM2", 9, exon),
           exon = ifelse(isotype == "CD", exon, gsub("8", "5", exon)),
           exon = ifelse(isotype == "CD", exon, gsub("9", "6", exon)),
           exon = gsub("[AB]", "", exon),
           exon = ifelse(strand == "+", exon, StrRev(exon)))
  )