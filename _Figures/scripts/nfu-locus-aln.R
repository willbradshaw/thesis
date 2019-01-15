#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

rm(list = ls())

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")
source("aux/blast.R")
source("aux/ggplot2.R")

bac_spread_path <- "/home/will/Downloads/cdown/summary_refiltered_spread.tsv"
nfu_aln_path <- "../_Data/locus/nfu_aln_blastn.tab"

blast_fmt <- "6 qseqid sseqid pident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send evalue bitscore qlen slen"

locus_name <- "IGH locus"
summary_single <- list(
  "165M01" = TRUE,
  "208A08" = TRUE,
  "277J10" = TRUE,
  "276N03" = FALSE,
  "scf10901" = TRUE,
  "scf9157" = FALSE,
  "209K12" = TRUE
)
summary_multi <- list(
  "chr6" = c(TRUE, 3)
)
summary_single[[locus_name]] <- TRUE

filename <- "nfu-locus-aln"
plot_height <- 15
plot_ratio <- 1.5

#------------------------------------------------------------------------------
# Prepare BAC locus table
#------------------------------------------------------------------------------

bac_spread <- suppressMessages(read_tsv(bac_spread_path)) %>% # Import table
  select(-matches("TM[12]")) %>% # Discount TM exons
  group_by(sseqid, slen) %>%
  summarise(vh = sum(vh_1), jh = sum(jh_1),
            cm = sum(ighm_1 + ighm_2 + ighm_3 + ighm_4),
            cd = sum(ighd_1 + ighd_2 + ighd_3 + ighd_4 +
                 ighd_5 + ighd_6 + ighd_7))

#------------------------------------------------------------------------------
# Prepare BLAST table
#------------------------------------------------------------------------------
qseqid <- "IGH"

btab <- read_blast_table(nfu_aln_path, blast_fmt)

nfu_aln <- btab %>%
  filter(length > 150) %>%
  mutate(sstart = sstart + as.integer(sapply(str_split(sseqid, "_"), function(x) x[2])) - 1,
         send = send + as.integer(sapply(str_split(sseqid, "_"), function(x) x[2])) - 1,
         slen = as.integer(sapply(str_split(sseqid, "_"), function(x) x[4])),
         sseqid = sapply(str_split(sseqid, "_"), function(x) x[1]))

qlen <- nfu_aln %>% pull(qlen) %>% head(1)

nfu_aln_concise <- nfu_aln %>%
  select(-qseqid, -qlen, -lstart, -lend, -nseq, -sstrand) %>%
  select(-mismatch, -gapopen, -gaps, -evalue)

#------------------------------------------------------------------------------
# Summarise alignments
#------------------------------------------------------------------------------

# Add row for locus itself
nfu_aln_concise <- tibble(sseqid = locus_name, pident = 100, qcovhsp = 100, 
                          length = qlen, qstart = 1, qend = qlen, sstart = 1, 
                          send = qlen, slen = qlen, orientation = TRUE, 
                          bitscore = 1e6) %>% # Made-up bitscore
  bind_rows(nfu_aln_concise)

# Get single-summary sequences
bitscore_threshold <- 1500
# Create summary tables for each sequence
summ_list_single <- list()
for (scf in names(summary_single)){
  ori <- summary_single[[scf]]
  tab <- nfu_aln_concise %>% 
    filter(sseqid == scf, orientation == ori, bitscore >= bitscore_threshold) %>%
    arrange(qstart)
  n_start <- which(tab$qstart == min(tab$qstart))
  n_end <- which(tab$qend == max(tab$qend))
  summ_list_single[[scf]] <- tab %>%
    summarise(sseqid = scf, orientation = ori, slen = first(slen),
              qstart = qstart[n_start], qend = qend[n_end],
              sstart = sstart[n_start], send = send[n_end]) %>%
    mutate(qstart_ext = qstart - ifelse(orientation, sstart - 1, slen - sstart),
           qend_ext = qend + ifelse(orientation, slen - send, send - 1),
           sstart_ext = 1, send_ext = slen, m_aln = 1)
}

# Combine single-match sequences and filter alignment table
summ_tab <- bind_rows(summ_list_single)
nfu_aln_cut <- filter(nfu_aln_concise, !(sseqid %in% names(summ_list_single)))

# Process multi-mapping regions
summ_list_multi <- list()
for (scf in names(summary_multi)){
  ori <- as.logical(summary_multi[[scf]][1])
  n_aln <- summary_multi[[scf]][2]
  tab <- nfu_aln_concise %>% 
    filter(sseqid == scf, orientation == ori, bitscore >= bitscore_threshold) %>%
    arrange(qstart) %>% 
    mutate(gap = qstart - lag(qend), gap = ifelse(is.na(gap), 0, gap),
           gap_rank = dense_rank(desc(gap)), new_aln = gap_rank < n_aln,
           m_aln = cumsum(new_aln)+1) %>%
    group_by(m_aln) %>%
    mutate(n_start = which(qstart == min(qstart)), n_end = which(qend == max(qend)))
  summ_list_multi[[scf]] <- tab %>%
    summarise(sseqid = scf, orientation = ori, slen = first(slen),
              qstart = nth(qstart, first(n_start)), qend = nth(qend, first(n_end)),
              sstart = nth(sstart, first(n_start)), send = nth(send, first(n_end))) %>%
    ungroup() %>%
    mutate(qstart_ext = qstart, qend_ext = qend, 
           sstart_ext = sstart, send_ext = send) # Don't extend these ones for now
}

# Re-combine into single summary table and add categories
summ_tab <- bind_rows(summ_tab, bind_rows(summ_list_multi)) %>%
  mutate(type = ifelse(sseqid == locus_name, "locus",
                       ifelse(grepl("scf", sseqid), "scf",
                              ifelse(grepl("chr", sseqid), "chr", "BAC")))) %>%
  arrange(factor(type, levels=c("locus", "BAC", "scf", "chr")))

# Extend chromosome ends to ends of plot
row_chr_start <- which(summ_tab$sseqid == "chr6" & summ_tab$m_aln == 1)
row_chr_end <- which(summ_tab$sseqid == "chr6" & summ_tab$m_aln == summary_multi[["chr6"]][2])
summ_tab[row_chr_start,]$qstart_ext <- min(summ_tab$qstart_ext)
summ_tab[row_chr_start,]$sstart_ext <- summ_tab[row_chr_start,]$sstart - 
  summ_tab[row_chr_start,]$qstart_ext + summ_tab[row_chr_start,]$qstart
summ_tab[row_chr_end,]$qend_ext <- max(summ_tab$qend_ext)
summ_tab[row_chr_end,]$send_ext <- summ_tab[row_chr_end,]$send + 
  summ_tab[row_chr_end,]$qend_ext - summ_tab[row_chr_end,]$qend

# Factorise sequence names for plotting
sseqid_factor <- factor(summ_tab$sseqid, levels = unique(summ_tab$sseqid))
summ_tab$sseqid <- sseqid_factor
level_chr <- which(levels(sseqid_factor) == "chr6")

#------------------------------------------------------------------------------
# Make plot
#------------------------------------------------------------------------------

chr_axis_breaks <- c(summ_tab %>% filter(sseqid == "chr6") %>% pull(qstart) %>% min,
                     summ_tab %>% filter(sseqid == "chr6") %>% pull(qend) %>% max)
chr_axis_labels <- c(summ_tab %>% filter(sseqid == "chr6") %>% pull(sstart) %>% min,
                     summ_tab %>% filter(sseqid == "chr6") %>% pull(send) %>% max)

g <- ggplot(summ_tab) + 
  geom_segment(aes(x = qstart_ext, xend = qend_ext,
                   y = factor(sseqid, levels=unique(sseqid)),
                   yend = factor(sseqid, levels=unique(sseqid)),
                   colour = orientation), 
               size = 5, lineend = "butt", linejoin = "mitre") + 
  xlab("Locus co-ordinate") + ylab("Sequence") +
  geom_segment(aes(x=min(qstart_ext-2e4), xend=max(qend_ext+2e4),
                   y = level_chr, yend = level_chr), linetype = 3,
               size = 1, colour = gg_color_hue(2)[2]) +
  scale_x_continuous(breaks = c(0, qlen),
                     sec.axis = dup_axis(name = "Chromosome co-ordinate",
                                         breaks = chr_axis_breaks,
                                         labels = chr_axis_labels)) +
  scale_colour_discrete(name = "Orientation",
                        labels = c("3' → 5'", "5' → 3'")) +
  theme_minimal() + theme_base +
  theme(axis.title.x = element_text(margin = margin(t=3, b=3, unit="mm")),
        axis.title.x.top = element_text(margin = margin(b=4, unit="mm")))

#------------------------------------------------------------------------------
# Save plot
#------------------------------------------------------------------------------

savefig(g, filename, height = plot_height, ratio = plot_ratio)
