###############################################################################
## AUX FILE                                                                  ##
## Auxiliary functions for range extraction                                  ##
###############################################################################

# Source packages
source("aux/packages.R")

#------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
#------------------------------------------------------------------------------

twostrand_match <- function(query, subject){
  # Check for exact matches of a pattern and its reverse complement in a sequence
  match_fwd <- vmatchPattern(query, subject)
  if (length(match_fwd[[1]])>0) return(match_fwd[[1]] %>% as.tibble %>% mutate(strand="+"))
  match_rev <- vmatchPattern(reverseComplement(query), subject)
  if (length(match_rev[[1]])>0) return(match_rev[[1]] %>% as.tibble %>% mutate(strand="-"))
  return(tibble())
}

twostrand_match_multi <- function(query, subjects){
  # Perform twostrand_match against a multisubject StringSet
  lapply(seq_along(subjects), function(n)
    twostrand_match(query, subjects[n]) %>% 
      mutate(seqnames = names(subjects)[n])
  ) %>% bind_rows
}

twostrand_alignment <- function(query, subject, type){
  # Align a query to a subject both ways and return the better alignment
  aln_fwd <- pairwiseAlignment(query, subject, type=type)
  aln_rev <- pairwiseAlignment(reverseComplement(query), subject, type=type)
  if (aln_fwd@score >= aln_rev@score) return(list(aln=aln_fwd, strand = "+"))
  return(list(aln=aln_rev, strand = "-"))
}

twostrand_alignment_multi <- function(query, subjects, type){
  alns <- lapply(seq_along(subjects), function(n)
    twostrand_alignment(query, subjects[n], type = type))
  scores <- sapply(alns, function(a) a$aln@score)
  return(alns[[which(scores == max(scores))[1]]])
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
  aln@subject@range %>% as.tibble %>%
    mutate(strand = strand,
           seqnames = names(aln@subject@unaligned)[1],
           label = names(aln@pattern@unaligned)[1])
}

align_mask <- function(query_seqs, subject_seqs, max_matches = 1,
                       type = "global-local"){
  # Successively align query sequences to a subject, masking each alignment
  # range before aligning next query
  aln_list <- list()
  range_list <- list()
  names_rejected <- character(0)
  names_out <- names(query_seqs)
  # Sort query seqs in descending order of length
  query_seqs <- query_seqs[order(width(query_seqs), decreasing = TRUE)]
  for (n in seq(length(query_seqs))){
    query_seq <- query_seqs[n]
    # First check for exact matches (much faster)
    aln <- twostrand_match_multi(query_seq[[1]], subject_seqs)
    if (nrow(aln) > max_matches){
      names_rejected <- c(names_rejected, names(query_seq))
      names_out <- names_out[!names_out %in% names(query_seq)]
    } else if (nrow(aln) > 0){ # Leave aln_list entry blank, fill in others
      range_list[[n]] <- aln[1,] %>% mutate(label = names(query_seq)[1])
      aln_list[[n]] <- list()
      subject_seqs[[range_list[[n]]$seqnames]][range_list[[n]]$start:range_list[[n]]$end] <- "N"
    } else { # Do full pairwise alignment (much slower)
      aln <- twostrand_alignment_multi(query_seq, subject_seqs, type = type)
      aln_list[[n]] <- aln$aln
      range_list[[n]] <- aln_to_range_tbl(aln$aln, aln$strand)
      subject_seqs[[range_list[[n]]$seqnames]][range_list[[n]]$start:range_list[[n]]$end] <- "N"
    }
  }
  names(aln_list) <- names_out
  names(range_list) <- names_out
  return(list(aln = aln_list, ranges = range_list, 
              subject_masked = subject_seqs, rejected = names_rejected))
}

align_match_groups <- function(match_names, query_seqs, subject_seq,
                               type = "global-local"){
  # Group query sequences into groups of identical sequences, find
  # matches in a subject sequence, and assign among grouped sequences
  # based on order of occurrence
  match_names <- match_names[order(match(match_names, names(query_seqs)))]
  aln_list <- list()
  range_list <- list()
  names_rejected <- character(0)
  while(length(match_names) > 0){
    # Identify first matching group
    match_group <- match_names[query_seqs[match_names] == query_seqs[match_names[1]]]
    match_seqs <- query_seqs[match_group]
    # Find n best alignments
    match_aln <- align_mask(match_seqs, subject_seq, type = type,
                            max_matches = length(match_seqs))
    # Assign alignments to group members based on order
    aln_order <- order(sapply(match_aln$ranges, function(r) r$start))
    names(match_aln$aln) <- match_group[aln_order]
    names(match_aln$ranges) <- match_group[aln_order]
    # Append alignments and ranges and move to next match group
    aln_list <- append(aln_list, match_aln$aln)
    range_list <- append(range_list, match_aln$ranges)
    subject_seq <- match_aln$subject_masked
    match_names <- match_names[! match_names %in% match_group]
    names_rejected <- c(names_rejected, match_aln$rejected)
  }
  return(list(aln = aln_list, ranges = range_list, 
              subject_masked = subject_seq, rejected = names_rejected))
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
  stop("Not all annotation entries processed into comments. ")
}