###############################################################################
## Generate CH range tables for multispecies locus maps                      ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

rm(list = ls())

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/io.R")
source("aux/ranges.R")

# Paths to genome sequences
species_ids <- c("ola", "cva", "fhe", "pfo", "pre", "kma", "ali", 
                 "ppl", "cto", "aau", "nor")
genome_paths <- lapply(species_ids, function(s) 
  paste0("../_Data/locus/scaffolds/", s, ".fasta"))
names(genome_paths) <- species_ids
genome_paths[["ola"]] <- sub("ola", "ola_new", genome_paths[["ola"]])

# Paths to constant exon sequences
segment_paths <- lapply(species_ids, function(s) 
  paste0("../_Data/constant/", s, "/", s, "_ch_nt.fasta"))
names(segment_paths) <- species_ids
segment_paths[["ola"]] <- sub("ola_ch", "ola_new_ch", segment_paths[["ola"]])

#------------------------------------------------------------------------------
# IMPORT SEQUENCES
#------------------------------------------------------------------------------

# Ordered segments
segments <- lapply(segment_paths, readDNAStringSet)

# Genomes
genomes <- lapply(genome_paths, readDNAStringSet)

# Strip terminal N's from segments
segments <- lapply(segments, function(t) 
  DNAStringSet(sapply(t, function(s) DNAString(gsub("^N+|N+$", "", s)))))

#------------------------------------------------------------------------------
# ALIGN SEGMENTS TO LOCUS AND EXTRACT RANGES
#------------------------------------------------------------------------------

# Identify pairs of identical sequences
names_matched <- lapply(segments, get_matched)
names_unmatched <- lapply(segments, get_unmatched)

# Align and mask singleton sequences
aln_unmatched <- list()
for (s in species_ids){ aln_unmatched[[s]] <- 
  align_mask(segments[[s]][names_unmatched[[s]]], genomes[[s]], 1)
}

# Align and mask matched sequences
aln_matched <- list()
for (s in species_ids){ aln_matched[[s]] <-
  align_match_groups(names_matched[[s]], segments[[s]],
                     aln_unmatched[[s]]$subject_masked) }

# Deal with rejected multimapping sequences (if any)
names_rejected <- lapply(species_ids, function(s)
  c(aln_matched[[s]]$rejected, aln_unmatched[[s]]$rejected))
names(names_rejected) <- species_ids
aln_rejected <- list()
for (s in species_ids){
  if (length(names_rejected[[s]]) == 0){
    aln_rejected[[s]] <- list(aln = list(), ranges = list(), 
                              rejected = list(), 
                              subject_masked = aln_matched[[s]]$subject_masked)
  } else {
    aln_rejected[[s]] <- align_mask(segments[[s]][names_rejected[[s]]], 
                                    aln_matched[[s]]$subject_masked, 1)
  }
}

n_still_rejected <- sapply(aln_rejected, function(t) length(t$rejected))

# Collate segment ranges
ranges <- lapply(names(segment_paths), function(s)
  append(append(aln_unmatched[[s]]$ranges, aln_matched[[s]]$ranges), 
         aln_rejected[[s]]$ranges) %>%
    bind_rows %>% 
    mutate(type = "CH",
           annotation = ifelse(!grepl("_", label), "", sub(".*?_(.*)", "\\1", label)),
           label = sub("_.*", "", label)) %>%
    arrange(seqnames, start, end)
)
names(ranges) <- names(segment_paths)
  
#------------------------------------------------------------------------------
# WRITE RANGES
#------------------------------------------------------------------------------

for (s in names(ranges)){
  d <- file.path("../_Data/ranges", s)
  if (!dir.exists(d)) dir.create(d)
  write_tsv(ranges[[s]] %>% select(-annotation),
            file.path(d, paste0(s, "-ch-ranges.tsv")))
}