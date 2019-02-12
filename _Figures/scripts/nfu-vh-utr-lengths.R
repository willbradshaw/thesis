###############################################################################
## TABLE                                                                     ##
## Get UTR sequence lengths from pilot IgSeq data                            ##
###############################################################################

tab_path <- "~/Downloads/cdown/seqs-all_mutations.tab" # TODO: Change this to somewhere in thesis dir
vh_name_conversion_path <- "../_Data/segments/nfu/nfu_vh_name_conversion.csv"

# Output path
outpath <- "../_Data/segments/nfu/nfu_vh_utr_lengths.tsv"

# Mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

# Read table
tab <- suppressMessages(read_tsv(tab_path))
vh_name_conversion <- suppressMessages(read_csv(vh_name_conversion_path))

# List possible V-assignments
vh_all <- tab %>% pull(BEST_V_CALL) %>% unique

# Exclude sequences with misplaced V-alignments and add UTR lengths
tab <- tab %>% filter(V_GERM_START_VDJ == 1) %>%
  mutate(UTR_LEN = V_SEQ_START - 1)

#------------------------------------------------------------------------------
# GET UTR LENGTH FOR CONSISTENT SEGMENTS
#------------------------------------------------------------------------------

# Get UTR length statistics
vh_utr_tab <- tab %>% group_by(BEST_V_CALL) %>%
  summarise(MEAN = mean(UTR_LEN),
            MEDIAN = median(UTR_LEN),
            MODE = Mode(UTR_LEN),
            SD = sd(UTR_LEN))

# Identify V-sequences where UTR mode matches median
vh_utr_tab_confident <- vh_utr_tab %>% filter(MEDIAN == MODE) %>%
  select(V = BEST_V_CALL, UTR_LEN = MEDIAN)
vh_confident <- vh_utr_tab_confident %>% pull(V)

#------------------------------------------------------------------------------
# ESTIMATE LENGTH FOR OTHER SEGMENTS
#------------------------------------------------------------------------------

# Filter input table to exclude confident Vs
tab <- tab %>% filter(!BEST_V_CALL %in% vh_confident)

window_size = 9
vh_utr_tab_unconfident <- tab %>% group_by(BEST_V_CALL) %>%
  mutate(UTR_RANK = row_number(UTR_LEN)) %>%
  arrange(UTR_RANK) %>%
  summarise(MODE = Mode(rollapply(UTR_LEN, width = window_size, 
                                  FUN = function(x) median(x, na.rm = TRUE)))) %>%
  filter(!is.na(MODE)) %>% select(V = BEST_V_CALL, UTR_LEN = MODE)

#------------------------------------------------------------------------------
# ADD NA ENTRIES FOR MISSING V-SEGMENTS
#------------------------------------------------------------------------------

vh_assigned <- c(vh_utr_tab_unconfident$V, vh_utr_tab_confident$V)
vh_missing <- vh_all[!vh_all %in% vh_assigned]
vh_utr_tab <- bind_rows(vh_utr_tab_confident,
                        vh_utr_tab_unconfident,
                        tibble(V = vh_missing, UTR_LEN = NA)) %>%
  rename(V_NAME_NEW = V) %>% 
  inner_join(vh_name_conversion, by = "V_NAME_NEW") %>%
  select(V_NAME_OLD, V_NAME_NEW, UTR_LEN) %>%
  arrange(V_NAME_OLD)

#------------------------------------------------------------------------------
# SAVE LENGTH ENTRIES
#------------------------------------------------------------------------------

write_tsv(vh_utr_tab, outpath)