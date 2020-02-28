#################################################################
## AUXILIARY R FUNCTIONS FOR PROCESSING BLAST ALIGNMENT TABLES ##
#################################################################

#-----------------------------------------------------------------------------
# Parse and process alignment data from a BLAST output table
#-----------------------------------------------------------------------------

split_format_string <- function(format){
    # Process a BLAST format string into a colnames vector
    fmt <- str_split(format, " ")[[1]]
    if (fmt[1] %in% as.character(seq(0,14))) fmt <- fmt[-1]
    return(fmt)
}

process_aln_table <- function(tab, extend_start = TRUE){
  # Process alignment co-ordinate columns from an alignment hit table
  # Required columns: orientation, sstart, send
  col_list <- c("orientation", "qstart", "qend", "sstart", "send")
  if (!all(c("orientation", "sstart", "send") %in% colnames(tab))){
    cols <- paste(col_list, collapse=", ")
    stop(paste("Required columns for table processing:", cols))
  }
  tab <- tab %>% mutate(lstart = sstart, lend = send, xstart = lstart, nseq=1)
  if (extend_start) tab <- tab %>% mutate(xstart = sstart - qstart + 1)
  tab[!tab$orientation,] <- tab[!tab$orientation,] %>%
    mutate(lstart = send, lend = 2*sstart - xstart)
  return(tab %>% select(-xstart))
}

read_blast_table <- function(path, format, delim = "\t"){
    # Read tabular blast output into a tibble
    tab <- suppressMessages(read_delim(path, delim=delim,
                      col_names = split_format_string(format)))
    if (nrow(tab) == 0) return(tab)
    return(tab %>%
             mutate(orientation = sstart < send) %>% 
             process_aln_table
    )
}

add_blast_regions <- function(blast_table, pattern, replacement){
    # Add region information for subsetting BLAST aligments
    if ("region" %in% names(blast_table)){
        blast_table %>% mutate(region = sub(pattern, replacement, region))
    } else {
        blast_table %>% mutate(region = sub(pattern, replacement, qseqid))
    }
}

intersect_region_ranges <- function(blast_table, regions){
    regions <- regions %>% transmute(region=label, rstart=start, rend=end)
    return(suppressMessages(inner_join(blast_table, regions)))
}

filter_by_regions <- function(blast_table, regions,
                             before, within, after, partial){
    # Filter BLAST alignments by the corresponding genomic ranges
    tab <- blast_table %>% intersect_region_ranges(regions)
    if (partial) {
        tab <- mutate(tab, tstart = lend, tend=lstart)
    } else {
        tab <- mutate(tab, tstart = lstart, tend = lend)
    }
    if (!before) tab <- filter(tab, tstart > rstart) # If tolerating overlaps, else lstart
    if (!within) tab <- filter(tab, !(tstart > rstart & tend < lend)) # If tol ol, else lstart/lend
    if (!after)  tab <- filter(tab, tend < rend) # If tol ol, else lend
    return(tab)
}

filter_by_quality <- function(blast_table){
    # Subset BLAST alignments by quality
    blast_table %>% group_by(qseqid) %>%
        filter(bitscore == max(bitscore)) %>% group_by()
}

filter_by_position <- function(blast_table, position, by_subject = FALSE){
    # Subset BLAST alignments by position
    if (position == "first"){
        f = min
    } else if (position == "last") {
        f = max
    } else {
        stop("Invalid position specifier: must be one of 'first' or 'last'.")
    }
    if (by_subject){
        blast_table <- group_by(blast_table, sseqid)
    } else {
        blast_table <- group_by(blast_table, qseqid)
    }
    return(blast_table %>% filter(lstart == f(lstart)) %>% ungroup())
}

filter_first <- function(blast_table, by_subject = FALSE){
    # Subset BLAST alignments to first in each matching set
    if (by_subject){
        blast_table <- group_by(blast_table, sseqid)
    } else {
        blast_table <- group_by(blast_table, qseqid)
    }
    return(blast_table %>% filter(row_number() == 1) %>% group_by())
}

blast_to_ranges <- function(blast_table, aa, trim_start=0, trim_end=0){
    # Convert a BLAST tibble into a GRanges object on a GenBank nt sequence
    if (nrow(blast_table) == 0) return(blast_table)
    return(blast_table %>% mutate(qlen = qlen * ifelse(aa, 3, 1)) %>%
           transmute(seqnames=sseqid, label=qseqid,
                     start=lstart + trim_start,
                     end = lstart + qlen - trim_end - 1,
                     width = end - start + 1,
                     strand = ifelse(orientation, "+", "-")))
}

add_types <- function(blast_table, type, subtype=""){
  # Add sequence type and subtype columns to an alignment table
  blast_table %>% mutate(type=type, subtype=subtype)
}

quality_threshold <- function(blast_table, e=NaN, b=NaN, i=NaN, n=NaN){
  if (!is.nan(e)) blast_table <- filter(blast_table, evalue <= e)
  if (!is.nan(b)) blast_table <- filter(blast_table, bitscore >= b)
  if (!is.nan(i)) blast_table <- filter(blast_table, pident >= i)
  if (!is.nan(n)) blast_table <- filter(blast_table, nseq >= n)
  return(blast_table)
}

group_overlaps <- function(blast_table, by_orientation = TRUE){
  # Group rows of BLAST table that overlap in their alignment ranges
  blast_table %>% arrange(lstart) %>%
    mutate(agroup=cumsum(cummax(lag(lend, default=dplyr::first(lend))) < lstart)) %>%
    group_by(agroup, add=TRUE) %>%
    {if(by_orientation) group_by(., orientation, add=TRUE) else .}
}

collapse_overlaps <- function(blast_table){
  # Collapse overlapping alignments in a BLAST table into a single sequence
  # range, retaining the best alignment score and counting sequences per group
  blast_table %>% group_overlaps() %>% 
    summarise(astart=min(lstart), aend=max(lend), nseq=sum(nseq),
              evalue = min(evalue, na.rm = TRUE),
              bitscore = max(bitscore, na.rm = TRUE)) %>%
    mutate(length = aend - astart) %>%
    ungroup() %>% arrange(sseqid, astart)
}

recollapse_overlaps <- function(overlap_table){
  # Re-collapse a table of collapsed overlap groups
  overlap_table %>% mutate(lstart=astart, lend=aend) %>%
    collapse_overlaps
}

regroup_overlaps <- function(overlap_table, by_orientation){
  # Re-group overlapping alignments from an overlap table
  overlap_table %>% mutate(lstart=astart, lend=aend) %>%
    group_overlaps(by_orientation = by_orientation) %>% 
    select(-lstart, -lend)
}
