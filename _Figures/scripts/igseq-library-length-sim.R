###############################################################################
## FIGURE                                                                    ##
## Expected library length distribution based on killifish V/D/J segments    ##
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
source("aux/ggplot2.R")

# Segment input paths
vh_nt_path <- "../_Data/segments/nfu/nfu_vh_nt_annotated.fasta"
dh_nt_path <- "../_Data/segments/nfu/nfu_dh_nt.fasta"
jh_nt_path <- "../_Data/segments/nfu/nfu_jh_nt.fasta"
ch_nt_path <- "../_Data/constant/nfu/nfu_ch_nt_exons.fasta"

# Path to VH UTR sequence lengths
vh_utr_lengths_path <- "../_Data/segments/nfu/nfu_vh_utr_lengths.tsv"

# Oligo input paths
primers_path <- "../_Data/oligos/igseq_primers.fasta"
tsa_path <- "../_Data/oligos/igseq_tsa.fasta"
p1_path <- "../_Data/oligos/truseq_adaptors_p1.fasta"
p2_path <- "../_Data/oligos/truseq_adaptors_p2.fasta"

# TODO: Paths to indel parameters here

#------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
#------------------------------------------------------------------------------

combine_vdjc <- function(vh, dh, jh, ch){
  # Concatenate all possible VDJC sequence combinations, retaining VH name
  vdh <- xscat(rep(vh, length(dh)), rep(dh, each = length(vh)))
  names(vdh) <- rep(names(vh), length(dh))
  vdjh <- xscat(rep(vdh, length(jh)), rep(jh, each = length(vdh)))
  names(vdjh) <- rep(names(vdh), length(jh))
  vdjch <- xscat(rep(vdjh, length(ch)), rep(ch, each = length(vdjh)))
  names(vdjch) <- rep(names(vdjh), length(ch))
  return(vdjch)
}

add_adaptors <- function(seqs, tsa, primers, p1, p2){
  # Add TSA and Illumina adaptor sequences to a DNAStringSet
  seqs_tsa <- xscat(tsa, seqs)
  m1s_ranges <- unlist(vmatchPattern(primers[["M1s"]], seqs_tsa))
  seqs_tsa_trimmed <- subseq(seqs_tsa, start(m1s_ranges))
  seqs_tsa_adapted <- xscat(p2[grepl("generic", names(p2))],
                            seqs_tsa_trimmed,
                            reverseComplement(p1[grepl("generic", names(p1))]))
  names(seqs_tsa_adapted) <- names(seqs)
  return(seqs_tsa_adapted)
}

ss_to_tibble <- function(seqs, vh_utr_lengths){
  # Convert stringSet to tibble and add UTR lengths
  tibble(V_NAME_OLD = names(seqs), WIDTH = width(seqs)) %>%
    inner_join(vh_utr_lengths, by = "V_NAME_OLD") %>%
    filter(!is.na(UTR_LEN)) %>%
    mutate(WIDTH = WIDTH + UTR_LEN) %>%
    select(-UTR_LEN)
}

#------------------------------------------------------------------------------
# IMPORT SEQUENCES
#------------------------------------------------------------------------------

# Import V/D/J segments
vh <- readDNAStringSet(vh_nt_path)
dh <- readDNAStringSet(dh_nt_path)
jh <- readDNAStringSet(jh_nt_path)
ch <- readDNAStringSet(ch_nt_path)

# Import UTR sequence length estimates
vh_utr_lengths <- suppressMessages(read_tsv(vh_utr_lengths_path))

# Import primers
tsa <- readDNAStringSet(tsa_path)
primers <- readDNAStringSet(primers_path)
p1 <- readDNAStringSet(p1_path)
p2 <- readDNAStringSet(p2_path)

# Filter CH sequences to just CM1 and truncate by primer sequence
ch <- ch[grepl("IGH.M-1", names(ch))]
cprimer_seq <- reverseComplement(primers[["CC"]])
primer_ranges <- unlist(vmatchPattern(cprimer_seq, ch))
ch <- subseq(ch, 1, end(primer_ranges))

#------------------------------------------------------------------------------
# PREPARE TABLE OF POSSIBLE SEQUENCE LENGTHS
#------------------------------------------------------------------------------

# Split segments by sublocus
vh_igh1 <- vh[grepl("IGH1", names(vh))]
dh_igh1 <- dh[grepl("IGH1", names(dh))]
jh_igh1 <- jh[grepl("IGH1", names(jh))]
ch_igh1 <- ch[grepl("IGH1", names(ch))]
vh_igh2 <- vh[grepl("IGH2", names(vh))]
dh_igh2 <- dh[grepl("IGH2", names(dh))]
jh_igh2 <- jh[grepl("IGH2", names(jh))]
ch_igh2 <- ch[grepl("IGH2", names(ch))]

# Generate all possible sequence combinations
vdjc <- c(combine_vdjc(vh_igh1, dh_igh1, jh_igh1, ch_igh1),
          combine_vdjc(vh_igh2, dh_igh2, jh_igh2, ch_igh2))

# Add TSA and adaptor sequences
vdjc_adapted <- add_adaptors(vdjc, tsa, primers, p1, p2)

# Convert to table, add UTR lengths and filter
vdjc_lengths <- ss_to_tibble(vdjc_adapted, vh_utr_lengths)

#------------------------------------------------------------------------------
# SUMMARISE BY V-IDENTITY AND GET RANGES
#------------------------------------------------------------------------------
vdjc_lengths_summary <- vdjc_lengths %>% 
  group_by(V_NAME_OLD) %>% 
  summarise(WMIN = min(WIDTH), WMAX = max(WIDTH))

g <- ggplot(vdjc_lengths_summary) + 
  geom_errorbar(aes(x = V_NAME_OLD, ymin = WMIN, ymax = WMAX)) +
  coord_flip() + xlab("VH segment") + ylab("Library sequence length (bp)") +
  theme_minimal() + theme_base


