###############################################################################
## FIGURE & DATA                                                             ##
## P20 and counts for V(D)J segment repertoires from pilot dataset           ##
###############################################################################
# TODO: Rename this
aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Output paths
filename_base <- "pilot-segments" # TODO: Fix this

# Set parameters
palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
             colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]])
n_exclude <- rep(5, 4)
n_exclude_null <- rep(0, 4)
names(n_exclude) <- paste0("2-0", seq(3,6))
names(n_exclude_null) <- paste0("2-0", seq(3,6))

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

tab <- import_tab(inpath_tab)
v_names <- suppressMessages(read_csv(inpath_vnames))
d_names <- suppressMessages(read_csv(inpath_dnames))
j_names <- suppressMessages(read_csv(inpath_jnames))

#------------------------------------------------------------------------------
# COUNT MISSING AND AMBIGUOUS SEGMENT CALLS
#------------------------------------------------------------------------------

# Count rates of missing and ambiguous segment calls
segtab_summary <- tab %>% 
  mutate(ONE_V = HAS_V & !is.na(V_AMBIG) & !V_AMBIG,
         ONE_D = HAS_D & !is.na(D_AMBIG) & !D_AMBIG,
         ONE_J = HAS_J & !is.na(J_AMBIG) & !J_AMBIG) %>%
  group_by(ONE_V, ONE_D, ONE_J) %>%
  summarise(N = n(), DUPCOUNT = sum(DUPCOUNT), CONSCOUNT = sum(CONSCOUNT)) %>%
  ungroup() %>%
  mutate(ONE_VJ = ONE_V & ONE_J, ONE_VDJ = ONE_VJ & ONE_D,
         N_PC = N/sum(N)*100, DUPCOUNT_PC = DUPCOUNT/sum(DUPCOUNT)*100,
         CONSCOUNT_PC = CONSCOUNT/sum(CONSCOUNT)*100)

# Save ambiguity measurements
savetxt(segtab_summary %>% filter(ONE_VJ) %>% pull(N_PC) %>% sum %>% round(1),
        outpath_one_vj_n_pc)
savetxt(segtab_summary %>% filter(ONE_VJ) %>% pull(CONSCOUNT_PC) %>% sum %>% round(1),
        outpath_one_vj_conscount_pc)
savetxt(segtab_summary %>% filter(ONE_VDJ) %>% pull(N_PC) %>% sum %>% round(1),
        outpath_one_vdj_n_pc)
savetxt(segtab_summary %>% filter(ONE_VDJ) %>% pull(CONSCOUNT_PC) %>% sum %>% round(1),
        outpath_one_vdj_conscount_pc)
savetxt(segtab_summary %>% filter(!ONE_D) %>% pull(N_PC) %>% sum %>% round(1),
        outpath_ambig_d_n_pc)
savetxt(segtab_summary %>% filter(!ONE_D) %>% pull(CONSCOUNT_PC) %>% sum %>% round(1),
        outpath_ambig_d_conscount_pc)

#------------------------------------------------------------------------------
# COMPUTE VJ RANKS AND SIZES
#------------------------------------------------------------------------------

vjtab <- tab %>% filter(!VJ_AMBIG) %>% 
  group_by(INDIVIDUAL, BEST_VJ_CALL, V_NAME_NEW = BEST_V_CALL, 
           J_NAME_NEW = BEST_J_CALL) %>%
  summarise(CLNCOUNT = n()) %>% group_by(INDIVIDUAL) %>% arrange(CLNCOUNT) %>%
  inner_join(v_names, by = "V_NAME_NEW") %>%
  inner_join(j_names, by = "J_NAME_NEW") %>%
  mutate(VJ_CALL = paste(V_NAME_OLD, J_NAME_OLD, sep = "/"),
         CLNRANK = row_number(desc(CLNCOUNT)),
         CLNFREQ = CLNCOUNT/sum(CLNCOUNT)) %>%
  arrange(CLNRANK) %>% ungroup %>%
  select(INDIVIDUAL, V_CALL = V_NAME_OLD, J_CALL = J_NAME_OLD,
         VJ_CALL, CLNCOUNT, CLNRANK, CLNFREQ)

# Rank:frequency plot
vj_rf <- zplot_line(vjtab, "INDIVIDUAL", palette, "Individual") +
  scale_x_log10(name = "VJ rank in repertoire") +
  scale_y_log10(name = "Relative VJ frequency")

# Save plot
ggsave(plot = vj_rf, filename = outpath_rf, device = "svg", units = "cm",
       height=15, width = 15*1.3)

#------------------------------------------------------------------------------
# COUNT VJ COMBINATIONS AND P20
#------------------------------------------------------------------------------

# P20
p20_vj <- vjtab %>% filter(CLNRANK <= 20) %>% group_by(INDIVIDUAL) %>% 
  summarise(P20 = sum(CLNFREQ) * 100)
savetxt(p20_vj %>% pull(P20) %>% min %>% round(1), "pilot-segments-vj-p20-min")
savetxt(p20_vj %>% pull(P20) %>% max %>% round(1), "pilot-segments-vj-p20-max")

# Theoretical and actual VJ combinations
n_vj_theoretical <- nrow(v_names) * nrow(j_names)
n_vj_any <- vjtab %>% group_by(VJ_CALL) %>% summarise %>%
  summarise(N = n()) %>% mutate(PC = N/n_vj_theoretical * 100)
n_vj <- vjtab %>% group_by(INDIVIDUAL, VJ_CALL) %>% summarise %>%
  summarise(N = n()) %>% mutate(PC = N/n_vj_theoretical * 100)

savetxt(n_vj_theoretical, outpath_n_vj_theoretical)
savetxt(n_vj_any %>% pull(N), outpath_n_vj_any)
savetxt(n_vj_any %>% pull(PC) %>% round(1), outpath_pc_vj_any)
savetxt(n_vj %>% pull(N) %>% min, outpath_n_vj_min)
savetxt(n_vj %>% pull(PC) %>% min %>% round(1), outpath_pc_vj_min)
savetxt(n_vj %>% pull(N) %>% max, outpath_n_vj_max)
savetxt(n_vj %>% pull(PC) %>% max %>% round(1), outpath_pc_vj_max)

#------------------------------------------------------------------------------
# SAVE P20 INFORMATION FOR COMPARISON WITH DIVERSITY SPECTRA
#------------------------------------------------------------------------------

# By individual
write_tsv(p20_vj, outpath_stats_indiv)

# By replicate
vjtab_rep <- tab %>% filter(!VJ_AMBIG) %>% 
  group_by(REPLICATE, BEST_VJ_CALL, V_NAME_NEW = BEST_V_CALL, 
           J_NAME_NEW = BEST_J_CALL) %>%
  summarise(CLNCOUNT = n()) %>% group_by(REPLICATE) %>% arrange(CLNCOUNT) %>%
  inner_join(v_names, by = "V_NAME_NEW") %>%
  inner_join(j_names, by = "J_NAME_NEW") %>%
  mutate(VJ_CALL = paste(V_NAME_OLD, J_NAME_OLD, sep = "/"),
         CLNRANK = row_number(desc(CLNCOUNT)),
         CLNFREQ = CLNCOUNT/sum(CLNCOUNT)) %>%
  arrange(CLNRANK) %>% ungroup %>%
  select(REPLICATE, V_CALL = V_NAME_OLD, J_CALL = J_NAME_OLD,
         VJ_CALL, CLNCOUNT, CLNRANK, CLNFREQ)
p20_vj_rep <- vjtab_rep %>% filter(CLNRANK <= 20) %>% group_by(REPLICATE) %>% 
  summarise(P20 = sum(CLNFREQ) * 100)
write_tsv(p20_vj_rep, outpath_stats_rep)
