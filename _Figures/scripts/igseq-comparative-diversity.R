###############################################################################
## FIGURE                                                                    ##
## Compare age-related decline in diversity between gut and body samples     ##
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
source("aux/changeo.R")

# Configure input
gut_clone_path <- "../_Data/changeo/spectra/gut-age_clone-diversity-solo_seqs-all_copy-NULL.tsv"
gut_vj_path <- "../_Data/changeo/spectra/gut-age_VJ-diversity-solo_seqs-all_copy-NULL.tsv"
age_clone_path <- "../_Data/changeo/spectra/ageing_clone-diversity-solo_seqs-all_copy-NULL.tsv"
age_vj_path <- "../_Data/changeo/spectra/ageing_VJ-diversity-solo_seqs-all_copy-NULL.tsv"
# TODO: Ageing data with large clones?

# Set parameters
qvals <- c(0,1,1.5,2,3,4)
significance_level <- 0.05

# Output paths
filename_base <- "igseq-comparative-diversity"

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

gut_clone <- import_div(gut_clone_path) %>% 
  mutate(SAMPLE = "Gut", TYPE = "Clonal repertoire")
gut_vj <- import_div(gut_vj_path) %>% 
  mutate(SAMPLE = "Gut", TYPE = "VJ repertoire")
age_clone <- import_div(age_clone_path) %>% 
  mutate(SAMPLE = "Whole body", TYPE = "Clonal repertoire")
age_vj <- import_div(age_vj_path) %>% 
  mutate(SAMPLE = "Whole body", TYPE = "VJ repertoire")

tab <- bind_rows(gut_clone, gut_vj, age_clone, age_vj) %>%
  mutate(AGE_DAYS = as.numeric(AGE_DAYS), AGE_WEEKS = as.numeric(AGE_WEEKS),
         AGE_DAYS = ifelse(is.na(AGE_DAYS), AGE_WEEKS*7, AGE_DAYS),
         AGE_WEEKS = ifelse(is.na(AGE_WEEKS), AGE_DAYS/7, AGE_WEEKS))


#------------------------------------------------------------------------------
# GENERATE AND COMPARE LINEAR MODELS / GLMs
#------------------------------------------------------------------------------


get_glm <- function(tab, q, type, family = Gamma(),
                    fl = formula("D~AGE_DAYS*SAMPLE")){
  tab %>% filter(Q == q, TYPE == type) %>%
  {if (is.na(family)) lm(fl, data = .)
    else glm(fl, family = family, data = .)}
}


tab_filtered <- tab %>% filter(Q == q)

tab_solo_filtered <- laggppply(qvals, function(q) tab_solo %>% filter(Q == q) %>%
                              select(AGE_DAYS, INDIVIDUAL, D) %>%

# Filter solo-diversity table
filter_solo_tab <- function(q){
  tab_solo %>% filter(Q == q) %>%
    select(AGE_DAYS, INDIVIDUAL, D) %>%
    mutate(AGE_DAYS = as.numeric(levels(AGE_DAYS))[AGE_DAYS]) %>%
    arrange(AGE_DAYS, INDIVIDUAL)
}
multi_filter <- function(qvals, family = Gamma()){
  bind_rows(lapply(qvals, function(q) 
    filter_solo_tab(q) %>% mutate(Q = q)))
}

# Generate (G)LMs
glm_diversity <- function(q, family = Gamma()){
  if (all(is.na(family))){
    lm(D~AGE_DAYS, data = filter_solo_tab(q))
  } else {
    glm(D~AGE_DAYS, family = family,
        data = filter_solo_tab(q))
  }
}
predict_glm_diversity <- function(q, family = Gamma()){
  g <- glm_diversity(q, family)
  ages <- tibble(AGE_DAYS = seq(min(as.numeric(age_groups)),
                                max(as.numeric(age_groups))))
  if (all(is.na(family))){
    p <- predict.lm(g, ages, type = "response")
  } else {
    p <- predict.glm(g, ages, type = "response")
  }
  return(bind_cols(ages, D = p))
}

multi_predict <- function(qvals, family = Gamma()){
  bind_rows(lapply(qvals, function(q) 
    predict_glm_diversity(q, family) %>% mutate(Q = q)))
}
