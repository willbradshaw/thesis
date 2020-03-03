###############################################################################
## ANALYSIS                                                                  ##
## Rarefaction of clonal counts for different datasets                       ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
#------------------------------------------------------------------------------

get_clonetab <- function(tab, scale, clone_field,
                         group_field){
  tab <- tab %>% filter(!is.na(!!as.name(clone_field)))
  if (is.null(scale)){
    clones <- tab %>% select(!!as.name(group_field), !!as.name(clone_field))
  } else {
    clones <- tibble(N = 1:sum(tab[[scale]]))
    clones[[group_field]] <- rep(tab[[group_field]], tab[[scale]])
    clones[[clone_field]] <- rep(tab[[clone_field]], tab[[scale]])
    clones <- select(clones, -N)
  }
  return(clones)
}

sample_clones_single <- function(tab, sample_size, scale,
                                 clone_field, group_field,
                                 replace, P){
  # Get individual-labelled vector of clone IDs
  clones <- get_clonetab(tab, scale, clone_field, group_field) %>%
    group_by(!!as.name(group_field)) %>%
    filter(n() >= sample_size)
  # Sample
  samples <- clones %>% sample_n(sample_size, replace = replace)
  # Count clones in new sample
  rarefied <- samples %>% 
    group_by(!!as.name(group_field), !!as.name(clone_field)) %>%
    summarise(N = n()) %>% 
    group_by(!!as.name(group_field)) %>%
    mutate(CLNRANK = row_number(desc(N)), CLNFREQ = N/sum(N),
           SINGLE = N == 1, SMALL = N < 5, LARGE = N >= 5)
  pcounts <- lapply(P, function(p)
    rarefied %>% filter(CLNRANK <= p) %>% 
      summarise(!!as.name(paste0("P", p)) := sum(CLNFREQ))) %>%
    join_all(by = "INDIVIDUAL", type = "full")
  rarefied_summ <- summarise(rarefied, N_CLONES = n(), 
                             N_CLONES_SINGLE = sum(SINGLE),
                             N_CLONES_SMALL = sum(SMALL),
                             N_CLONES_LARGE = sum(LARGE)) %>%
    mutate(PC_CLONES_SMALL = N_CLONES_SMALL/N_CLONES) %>%
    inner_join(pcounts, by = "INDIVIDUAL")
  return(rarefied_summ)
}

rarefy_clones_single <- function(tab, sample_size, n_repeats,
                                 scale, clone_field,
                                 group_field, replace, P){
  # Repeatedly downsample and obtain clone counts for a single
  # sample size
  rarefactions <- lapply(1:n_repeats, function(m)
    sample_clones_single(tab, sample_size, scale, 
                         clone_field, group_field, replace, P) %>%
    mutate(ITER = m)) %>%
    bind_rows %>% melt(id.vars = c("INDIVIDUAL", "ITER"), 
                       variable.name = "METRIC", value.name = "VALUE")
  rarefactions_summ <- rarefactions %>% 
    group_by(!!as.name(group_field), METRIC) %>%
    summarise(MEAN = mean(VALUE), SD = sd(VALUE)) %>%
    mutate(SAMPLE_SIZE = sample_size)
  return(rarefactions_summ)
}

rarefy_clones <- function(tab, sample_sizes, n_repeats, scale, P,
                          clone_field = "CLONE",
                          group_field = "INDIVIDUAL",
                          replace = FALSE){
  # Repeatedly downsample and obtain clone counts for a range
  # of sample sizes
  bind_rows(lapply(sample_sizes, function(s)
    rarefy_clones_single(tab, s, n_repeats, scale, clone_field, group_field,
                         replace, P)))
}

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

write_log("Importing data...", newline = FALSE)
tab_gut <- import_tab(inpath_gut) %>% 
  mutate(GROUP = sub("YI_16", "YI_6", GROUP)) %>%
  filter(!INDIVIDUAL %in% individuals_excluded)
tab_age <- import_tab(inpath_age)
tab_pilot <- import_tab(inpath_pilot)
log_done(newline=TRUE)

#------------------------------------------------------------------------------
# RAREFY CLONAL COUNTS
#------------------------------------------------------------------------------

write_log("Rarefying pilot dataset...", newline = FALSE)
r_pilot <- rarefy_clones(tab_pilot, sample_sizes, n_repeats, scale, p) %>%
  mutate(EXPERIMENT = "pilot")
log_done(newline=TRUE)
write_log("Rarefying ageing dataset...", newline = FALSE)
r_age <- rarefy_clones(tab_age, sample_sizes, n_repeats, scale, p) %>%
  mutate(EXPERIMENT = "ageing")
log_done(newline=TRUE)
write_log("Rarefying gut dataset...", newline = FALSE)
r_gut <- rarefy_clones(tab_gut, sample_sizes, n_repeats, scale, p) %>%
  mutate(EXPERIMENT = "gut")
log_done(newline=TRUE)

r <- bind_rows(r_gut, r_age, r_pilot)

#------------------------------------------------------------------------------
# SAVE OUTPUT
#------------------------------------------------------------------------------

write_tsv(r, outpath)
