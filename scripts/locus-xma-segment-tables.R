###############################################################################
## Xiphophorus maculatus new locus - VH range table inference                ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
#------------------------------------------------------------------------------

twostrand_match <- function(query, subject){
  # Check for exact matches of a pattern and its reverse complement in a sequence
  match_fwd <- vmatchPattern(query, subject)
  if (length(match_fwd[[1]])>0) return(match_fwd[[1]] %>% as_tibble %>% mutate(strand="+"))
  match_rev <- vmatchPattern(reverseComplement(query), subject)
  if (length(match_rev[[1]])>0) return(match_rev[[1]] %>% as_tibble %>% mutate(strand="-"))
  return(tibble())
}

twostrand_alignment <- function(query, subject, type = "global-local"){
  # Align a query to a subject both ways and return the better alignment
  aln_fwd <- pairwiseAlignment(query, subject, type=type)
  aln_rev <- pairwiseAlignment(reverseComplement(query), subject, type=type)
  if (aln_fwd@score >= aln_rev@score) return(list(aln=aln_fwd, strand = "+"))
  return(list(aln=aln_rev, strand = "-"))
}

get_matched <- function(stringset){
  # Get a vector of names of sequences from a StringSet that exactly match 
  # another member
  seg_match <- outer(stringset, stringset, `==`) %>%
    melt(varnames = c("seg1", "seg2"), value.name = "match") %>%
    filter(seg1 != seg2, match) %>%
    mutate(seg1 = as.character(seg1), seg2 = as.character(seg2))
  n <- 1
  while (n < nrow(seg_match)){
    seg_match <- seg_match %>% filter(seg1 != nth(seg2, n) | seg2 != nth(seg1, n))
    n <- n + 1
  }
  match_names <- unique(c(seg_match$seg1, seg_match$seg2))
  match_names <- match_names[order(match(match_names, names(stringset)))]
  return(match_names)
}

get_unmatched <- function(stringset){
  # Get a vector of names of sequences from a StringSet that do not exactly 
  # match any other member
  match_names <- get_matched(stringset)
  return(names(stringset)[! names(stringset) %in% match_names])
}

aln_to_range_tbl <- function(aln, strand){
  # Convert an alignment (plus strand info) to a range tibble
  aln@subject@range %>% as_tibble %>%
    mutate(strand = strand,
           seqnames = names(aln@subject@unaligned)[1],
           label = names(aln@pattern@unaligned)[1])
}

align_mask <- function(query_seqs, subject_seq, type = "global-local"){
  # Successively align query sequences to a subject, masking each alignment
  # range before aligning next query
  aln_list <- list()
  range_list <- list()
  # Sort query seqs in descending order of length
  query_seqs <- query_seqs[order(width(query_seqs), decreasing = TRUE)]
  for (n in seq(length(query_seqs))){
    query_seq <- query_seqs[n]
    # First check for exact matches (much faster)
    aln <- twostrand_match(query_seq[[1]], subject_seq)
    if (nrow(aln) > 0){ # Leave aln_list entry blank, fill in others
      range_list[[n]] <- aln[1,] %>% 
        mutate(seqnames = names(subject_seq)[1], label = names(query_seq)[1])
      aln_list[[n]] <- list()
    } else { # Do full pairwise alignment (much slower)
      aln <- twostrand_alignment(query_seq, subject_seq, type = type)
      aln_list[[n]] <- aln$aln
      range_list[[n]] <- aln_to_range_tbl(aln$aln, aln$strand)
    }
    subject_seq[[1]][range_list[[n]]$start:range_list[[n]]$end] <- "N"
  }
  names(aln_list) <- names(query_seqs)
  names(range_list) <- names(query_seqs)
  return(list(aln = aln_list, ranges = range_list, subject_masked = subject_seq))
}

align_match_groups <- function(match_names, query_seqs, subject_seq,
                               type = "global-local"){
  # Group query sequences into groups of identical sequences, find
  # matches in a subject sequence, and assign among grouped sequences
  # based on order of occurrence
  match_names <- match_names[order(match(match_names, names(query_seqs)))]
  aln_list <- list()
  range_list <- list()
  while(length(match_names) > 0){
    # Identify first matching group
    match_group <- match_names[query_seqs[match_names] == query_seqs[match_names[1]]]
    match_seqs <- query_seqs[match_group]
    # Find n best alignments
    match_aln <- align_mask(match_seqs, subject_seq, type = type)
    # Assign alignments to group members based on order
    aln_order <- order(sapply(match_aln$ranges, function(r) r$start))
    names(match_aln$aln) <- match_group[aln_order]
    names(match_aln$ranges) <- match_group[aln_order]
    # Append alignments and ranges and move to next match group
    aln_list <- append(aln_list, match_aln$aln)
    range_list <- append(range_list, match_aln$ranges)
    subject_seq <- match_aln$subject_masked
    match_names <- match_names[! match_names %in% match_group]
  }
  return(list(aln = aln_list, ranges = range_list, subject_masked = subject_seq))
}

rss_to_tbl <- function(rss_stringset, orientation = TRUE, dh = FALSE){
  # Convert a StringSet of RSS sequences to an output tibble
  rss_stringset <- rss_stringset[width(rss_stringset) > 0]
  heptamer_start <- ifelse(orientation, 1, -7)
  heptamer_end <- ifelse(orientation, 7, -1)
  nonamer_start <- ifelse(orientation, -9, 1)
  nonamer_end <- ifelse(orientation, -1, 9)
  spacer_start <- ifelse(orientation, heptamer_end + 1, nonamer_end + 1)
  spacer_end <- ifelse(orientation, nonamer_start - 1, heptamer_start - 1)
  heptamers <- subseq(rss_stringset, heptamer_start, heptamer_end)
  nonamers <- subseq(rss_stringset, nonamer_start, nonamer_end)
  spacers <- subseq(rss_stringset, spacer_start, spacer_end)
  if(dh) rss_type = sub(".*_rss_", "", names(rss_stringset))
  rss_tbl <- tibble(label = sub("_.*", "", names(rss_stringset)),
                    heptamer = as.character(heptamers),
                    spacer = as.character(spacers),
                    spacer_length = width(spacers),
                    nonamer = as.character(nonamers)
  )
  if (dh) rss_tbl <- mutate(rss_tbl, rss_type = rss_type)
  return(rss_tbl)
}

char_rc <- function(char){
  # Reverse-complement a character vector, keeping the format
  char %>% DNAStringSet %>% reverseComplement %>% as.character
}

purge_pattern <- function(pattern, vector, sep="_"){
  # Strip a pattern from a character vector, plus any flanking separators
  vector %>% gsub(pattern, "", .) %>% gsub(paste0(sep, sep), sep, .) %>%
    sub(paste0("^", sep), "", .) %>% sub(paste0(sep, "$"), "", .)
}

capitalise <- function(string){
  # Naively capitalise first character of each string in a character vector
  paste0(toupper(substring(string, 1, 1)), substring(string, 2))
}


annot_to_comment <- function(table, pattern, new_comment){
  # Extract annotations into formatted comments
  table %>%
    mutate(comment = ifelse(grepl(pattern, annotation),
                            ifelse(comment == "", capitalise(new_comment),
                                   paste(comment, new_comment, sep=", ")),
                            comment),
           annotation = purge_pattern(pattern, annotation)
    )
}

check_annotations <- function(table){
  # Confirm that annotation column has been fully converted to comments
  annots_remaining <- table$annotation[table$annotation != ""]
  if (length(annots_remaining) == 0) return(table %>% select(-annotation))
  stop(paste0("Not all annotation entries processed into comments. ",
              "Entries remaining: ",
              a, "."))
}

#------------------------------------------------------------------------------
# IMPORT SEQUENCES
#------------------------------------------------------------------------------

# Ordered segments
vh_nt <- readDNAStringSet(inpath_vh_nt)
dh_nt <- readDNAStringSet(inpath_dh_nt)
jh_nt <- readDNAStringSet(inpath_jh_nt)
jh_aa <- readAAStringSet(inpath_jh_aa)
ch_nt <- readDNAStringSet(inpath_ch_nt)

# RSS
vh_rss <- readDNAStringSet(inpath_vh_rss)
dh_rss <- readDNAStringSet(inpath_dh_rss)
jh_rss <- readDNAStringSet(inpath_jh_rss)

# Locus
locus_unmasked <- readDNAStringSet(inpath_locus)

#------------------------------------------------------------------------------
# ALIGN SEGMENTS TO LOCUS AND EXTRACT RANGES
#------------------------------------------------------------------------------

# Combine DH segments with RSS sequences for more unique sequence
dh_rss_tab <- tibble(label = names(dh_rss), sequence = as.character(dh_rss)) %>%
  mutate(type = sub(".*_rss_", "", label), seqnames = sub("_rss_.*", "", label)) %>%
  group_by(seqnames) %>% 
  summarise(sequence_3prime = sequence[type == "3prime"],
            sequence_5prime = sequence[type == "5prime"],
            sequence_5prime = as.character(reverseComplement(DNAString(sequence_5prime)))
            ) %>%
  inner_join(tibble(seqnames = names(dh_nt), sequence = as.character(dh_nt)),
                    by = "seqnames") %>%
  mutate(sequence_all = paste0(sequence_5prime, sequence, sequence_3prime))
dh_plus_rss <- DNAStringSet(dh_rss_tab$sequence_all)
names(dh_plus_rss) <- dh_rss_tab$seqnames
  
# Identify pairs of identical sequences
vh_matched <- get_matched(vh_nt)
vh_unmatched <- get_unmatched(vh_nt)
dh_matched <- get_matched(dh_plus_rss)
dh_unmatched <- get_unmatched(dh_plus_rss)
jh_matched <- get_matched(jh_nt)
jh_unmatched <- get_unmatched(jh_nt)
ch_matched <- get_matched(ch_nt)
ch_unmatched <- get_unmatched(ch_nt)

# Align and mask singleton sequences
vh_aln_unmatched <- align_mask(vh_nt[vh_unmatched], locus_unmasked)
dh_aln_unmatched <- align_mask(dh_plus_rss[dh_unmatched], locus_unmasked)
jh_aln_unmatched <- align_mask(jh_nt[jh_unmatched], locus_unmasked)
ch_aln_unmatched <- align_mask(ch_nt[ch_unmatched], locus_unmasked)

# Align and mask matched sequences
vh_aln_matched <- align_match_groups(vh_matched, vh_nt, vh_aln_unmatched$subject_masked)
jh_aln_matched <- align_match_groups(jh_matched, jh_nt, dh_aln_unmatched$subject_masked)
dh_aln_matched <- align_match_groups(dh_matched, dh_nt, jh_aln_unmatched$subject_masked)
ch_aln_matched <- align_match_groups(ch_matched, ch_nt, ch_aln_unmatched$subject_masked)

# Collate segment ranges
vh_ranges <- append(vh_aln_unmatched$ranges, vh_aln_matched$ranges) %>%
  bind_rows() %>% 
  mutate(type = "VH",
         annotation = ifelse(!grepl("_", label), "", sub(".*?_(.*)", "\\1", label)),
         label = sub("_.*", "", label)) %>%
  arrange(start)
dh_ranges <- append(dh_aln_unmatched$ranges, dh_aln_matched$ranges) %>%
  bind_rows() %>% 
  inner_join(dh_rss_tab %>% rename(label = seqnames), by = "label") %>%
  mutate(start = start + str_count(sequence_5prime),
         end = end - str_count(sequence_3prime),
         width = end - start + 1) %>% # Remove RSSs from ranges
  select(-contains("sequence")) %>%
  mutate(type = "DH",
         annotation = ifelse(!grepl("_", label), "", sub(".*?_(.*)", "\\1", label)),
         label = sub("_.*", "", label)) %>%
  arrange(start)
jh_ranges <- append(jh_aln_unmatched$ranges, jh_aln_matched$ranges) %>%
  bind_rows() %>% 
  mutate(type = "JH",
         annotation = ifelse(!grepl("_", label), "", sub(".*?_(.*)", "\\1", label)),
         label = sub("_.*", "", label)) %>%
  arrange(start)
ch_ranges <- append(ch_aln_unmatched$ranges, ch_aln_matched$ranges) %>%
  bind_rows() %>% 
  mutate(type = "CH",
         annotation = ifelse(!grepl("_", label), "", sub(".*?_(.*)", "\\1", label)),
         label = sub("_.*", "", label)) %>%
  arrange(start)

# TODO: Make printable table
# TODO: Repeat with other segment types

#------------------------------------------------------------------------------
# EXTRACT RSS INFORMATION
#------------------------------------------------------------------------------

vh_rss_tbl <- rss_to_tbl(vh_rss)
dh_rss_tbl <- rss_to_tbl(dh_rss, dh = TRUE)
jh_rss_tbl <- rss_to_tbl(jh_rss)

dh_rss5_tbl <- filter(dh_rss_tbl, rss_type == "5prime") %>% select(-rss_type)
dh_rss3_tbl <- filter(dh_rss_tbl, rss_type == "3prime") %>% select(-rss_type)

#------------------------------------------------------------------------------
# GENERATE SEGMENT TABLES FROM RANGES
#------------------------------------------------------------------------------

# Join segment and RSS information
vh_ranges_plus <- full_join(vh_ranges, vh_rss_tbl, by = "label") %>%
  mutate(rss_start = end + 1, 
         rss_width = str_count(heptamer) + str_count(nonamer) + spacer_length,
         rss_end = rss_start + rss_width - 1,
         comment = "")

jh_ranges_plus <- full_join(jh_ranges, jh_rss_tbl, by = "label") %>%
  mutate(rss_end = start - 1, 
         rss_width = str_count(heptamer) + str_count(nonamer) + spacer_length,
         rss_start = end - rss_width + 1,
         comment = "",
         heptamer = char_rc(heptamer),
         nonamer = char_rc(nonamer)
         )

dh_ranges_plus <- full_join(dh_ranges, dh_rss3_tbl, by = "label") %>%
  mutate(rss3_start = end + 1, 
         rss3_width = str_count(heptamer) + str_count(nonamer) + spacer_length,
         rss3_end = rss3_start + rss3_width - 1,
         comment = "") %>%
  rename(rss3_heptamer = heptamer, rss3_nonamer = nonamer, rss3_spacer = spacer,
         rss3_spacer_length = spacer_length) %>%
  full_join(dh_rss5_tbl, by = "label") %>%
  mutate(rss5_end = start - 1, 
         rss5_width = str_count(heptamer) + str_count(nonamer) + spacer_length,
         rss5_start = rss5_end - rss5_width + 1,
         heptamer = char_rc(heptamer),
         nonamer = char_rc(nonamer)) %>%
  rename(rss5_heptamer = heptamer, rss5_nonamer = nonamer, rss5_spacer = spacer,
         rss5_spacer_length = spacer_length)

ch_pattern <- "IGH.?([MZD]).?-(S?(TM)?\\d?)(.*)"
ch_ranges_plus <- ch_ranges %>%
  mutate(isotype = sub(ch_pattern, "\\1", label),
         exon = sub(ch_pattern, "\\2", label),
         variant = sub(ch_pattern, "\\4", label),
         comment = "")

# Parse annotations
annots <- list(
  "frameshift" = "frameshift",
  "nonsense" = "nonsense mutation",
  "truncated_3prime" = "3'-truncated",
  "truncated_5prime" = "5'-truncated",
  "no_rss" = "no RSS",
  "new" = "named from new locus (rename before publication)",
  "old" = "named from old locus (rename before publication)"
)
for (a in names(annots)){
  vh_ranges_plus <- annot_to_comment(vh_ranges_plus, a, annots[[a]])
  jh_ranges_plus <- annot_to_comment(jh_ranges_plus, a, annots[[a]])
  dh_ranges_plus <- annot_to_comment(dh_ranges_plus, a, annots[[a]])
  ch_ranges_plus <- annot_to_comment(ch_ranges_plus, a, annots[[a]])
}
vh_ranges_plus <- check_annotations(vh_ranges_plus)
dh_ranges_plus <- check_annotations(dh_ranges_plus)
jh_ranges_plus <- check_annotations(jh_ranges_plus)
ch_ranges_plus <- check_annotations(ch_ranges_plus)

# Add sequences to DH and JH tables
dh_ranges_plus <- dh_ranges_plus %>% 
  inner_join(tibble(label = names(dh_nt) %>% sub("_.*", "", .),
                    sequence_nt = as.character(dh_nt)), by = "label")
jh_ranges_plus <- jh_ranges_plus %>% 
  inner_join(tibble(label = names(jh_nt) %>% sub("_.*", "", .),
                    sequence_nt = as.character(jh_nt)), by = "label") %>%
  inner_join(tibble(label = names(jh_aa) %>% sub("_.*", "", .),
                    sequence_aa = as.character(jh_aa)), by = "label")

# Arrange column & add readable column names
vh_tab <- vh_ranges_plus %>% 
  mutate(heptamer = ifelse(is.na(heptamer), "", heptamer),
         nonamer = ifelse(is.na(nonamer), "", nonamer),
         spacer = ifelse(is.na(spacer), "", spacer),
         spacer_length = ifelse(is.na(spacer_length), "", spacer_length)) %>%
  select("Name" = label, "Start" = start, "End" = end, "Length" = width,
         "Strand" = strand,
         "RSS Start" = rss_start, "Heptamer" = heptamer, 
         "Spacer Length" = spacer_length, "Nonamer" = nonamer, 
         "RSS End" = rss_end, "RSS Length" = rss_width, "Comment" = comment)

jh_tab_seg <- jh_ranges_plus %>%
  select("Name" = label, "Start" = start, "NT Sequence" = sequence_nt,
         "AA Sequence" = sequence_aa, "End" = end, "Length" = width,
         "Strand" = strand, "Comment" = comment)
if (all(jh_tab_seg["Comment"] == "")) jh_tab_seg = select(jh_tab_seg, -`Comment`)

jh_tab_rss <- jh_ranges_plus %>% 
  mutate(heptamer = ifelse(is.na(heptamer), "", heptamer),
         nonamer = ifelse(is.na(nonamer), "", nonamer),
         spacer = ifelse(is.na(spacer), "", spacer),
         spacer_length = ifelse(is.na(spacer_length), "", spacer_length)) %>%
  select("Name" = label, "RSS Start" = rss_start, "Nonamer" = nonamer, 
         "Spacer Length" = spacer_length, "Heptamer" = heptamer,
         "RSS End" = rss_end, "RSS Length" = rss_width)

dh_tab_seg <- dh_ranges_plus %>%
  select("Name" = label, "Start" = start, "NT Sequence" = sequence_nt,
         "End" = end, "Length" = width,
         "Strand" = strand, "Comment" = comment)
if (all(dh_tab_seg["Comment"] == "")) dh_tab_seg = select(dh_tab_seg, -`Comment`)

dh_tab_rss5 <- dh_ranges_plus %>% 
  mutate(rss5_heptamer = ifelse(is.na(rss5_heptamer), "", rss5_heptamer),
         rss5_nonamer = ifelse(is.na(rss5_nonamer), "", rss5_nonamer),
         rss5_spacer = ifelse(is.na(rss5_spacer), "", rss5_spacer),
         rss5_spacer_length = ifelse(is.na(rss5_spacer_length), "", rss5_spacer_length)
         ) %>%
  select("Name" = label, "5'-RSS Start" = rss5_start, "Nonamer" = rss5_nonamer, 
         "Spacer Length" = rss5_spacer_length, "Heptamer" = rss5_heptamer,
         "5'-RSS End" = rss5_end, "Length" = rss5_width)

dh_tab_rss3 <- dh_ranges_plus %>% 
  mutate(rss3_heptamer = ifelse(is.na(rss3_heptamer), "", rss3_heptamer),
         rss3_nonamer = ifelse(is.na(rss3_nonamer), "", rss3_nonamer),
         rss3_spacer = ifelse(is.na(rss3_spacer), "", rss3_spacer),
         rss3_spacer_length = ifelse(is.na(rss3_spacer_length), "", rss3_spacer_length)) %>%
  select("Name" = label, "3'-RSS Start" = rss3_start, "Heptamer" = rss3_heptamer,
         "Spacer Length" = rss3_spacer_length, "Nonamer" = rss3_nonamer, 
         "3'-RSS End" = rss3_end, "Length" = rss3_width)

ch_tab <- ch_ranges_plus %>% 
  select("Name" = label, "Isotype" = isotype, "Start" = start, "End" = end, 
         "Length" = width, "Strand" = strand, "Comment" = comment)
if(all(ch_tab$Comment == "")) ch_tab <- select(ch_tab, -Comment)

#------------------------------------------------------------------------------
# WRITE RANGES
#------------------------------------------------------------------------------

#rangedir <- "../_Data/ranges/xma"
#for (segtype in paste0(c("v", "d", "j", "c"), "h")){
#  write_tsv(get(paste0(segtype, "_ranges")) %>% select(-annotation),
#            file.path(rangedir, get(paste0("rangename_", segtype))))
#}
#write_tsv(bind_rows(vh_ranges, dh_ranges, jh_ranges, ch_ranges) %>% 
#            select(-annotation), file.path(rangedir, rangename_all))

#------------------------------------------------------------------------------
# WRITE TABLES
#------------------------------------------------------------------------------

# Split VH table to fit on page
vh_per_tab <- 25
vh_split_at <- seq(0, ceiling(nrow(vh_tab)/vh_per_tab)*vh_per_tab, vh_per_tab)
vh_n_tab <- length(vh_split_at) - 1

for (n in 1:vh_n_tab){
  assign(paste0("vh_tab_", n), vh_tab[(vh_split_at[n]+1):vh_split_at[n+1],])
  savetab(get(paste0("vh_tab_", n)), get(paste0("outpath_vh", n)))
}

# Save to latex tables
savetab(dh_tab_seg, outpath_dh_seg)
savetab(dh_tab_rss5, outpath_dh_rss5)
savetab(dh_tab_rss3, outpath_dh_rss3)
savetab(jh_tab_seg, outpath_jh_seg)
savetab(jh_tab_rss, outpath_jh_rss)
savetab(ch_tab, outpath_ch)

#------------------------------------------------------------------------------
# DETERMINE % OF RSS WITH SPACER WITHIN 1BP OF EXPECTED
#------------------------------------------------------------------------------

#rss_match_vh <- vh_tab %>% 
#  mutate(spacer_diff = abs(as.numeric(`Spacer Length`) - 23)) %>%
#  pull(spacer_diff) %>% .[!is.na(.)]
#rss_match_jh <- jh_tab_rss %>% 
#  mutate(spacer_diff = abs(as.numeric(`Spacer Length`) - 23)) %>%
#  pull(spacer_diff) %>% .[!is.na(.)]
#rss_match_dh <- bind_rows(dh_tab_rss5, dh_tab_rss3) %>%
#  mutate(spacer_diff = abs(as.numeric(`Spacer Length`) - 12)) %>%
#  pull(spacer_diff) %>% .[!is.na(.)]

#rss_match_pc <- mean(c(rss_match_vh, rss_match_dh, rss_match_jh) <= 1) * 100
