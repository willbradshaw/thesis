###############################################################################
## TABLE                                                                     ##
## 2017 gut-microbiota study cohort fish and group summary                   ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

rm(list = ls())

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/io.R")

# Path to input data
inpath <- "../_Data/fish_data/gut-samples.csv"

# Output paths
outfile_fish <- "gut-cohorts-fish"
outfile_summary <- "gut-cohorts-summary"

#------------------------------------------------------------------------------
# IMPORT DATA AND GENERATE SUMMARY
#------------------------------------------------------------------------------

tab_in <- suppressMessages(read_csv(inpath))
# TODO: Add hatch and death dates, age in days, death weight data if available

treatments <- list(
  "SMT" = "Antibiotics + microbiota transfer from middle-aged fish",
  "YMT" = "Antibiotics + microbiota transfer from young fish",
  "ABX" = "Antibiotics treatment only",
  "YI" = "Young fish, no treatment",
  "WT" = "Old fish, no treatment"
)
antibiotics <- list(
  "SMT" = "Yes",
  "YMT" = "Yes",
  "ABX" = "Yes",
  "YI" = "No",
  "WT" = "No"
)
transfer <- list(
  "SMT" = "Yes (9.5-week-old donor)",
  "YMT" = "Yes (6-week-old donor)",
  "ABX" = "No",
  "YI" = "No",
  "WT" = "No"
)

sort_levels <- c("YI", "WT", "ABX", "SMT", "YMT")

tab_summary <- tab_in %>% group_by(GROUP, AGE_WEEKS, TREATMENT) %>%
  summarise(`# Fish (Sequenced/Total)` = paste(sum(LP > 0), n(), sep = " / "),
            #`Hatch date` = as.character(first(`Hatch date`)),
            #`Sacrifice date` = as.character(first(`Sacrifice date`)),
            #`Age (days)` = format(first(`Age (days)`), digits = 0),
            #`Age (weeks)` = format(first(`Age (weeks)`), digits = 3),
            #`Mean weight (g)` = format(mean(`Death weight (g)`), digits = 3)
            ) %>% ungroup %>%
  mutate(`Antibiotics?` = unlist(antibiotics[TREATMENT]),
         `Microbiota Transfer?` = unlist(transfer[TREATMENT])) %>%
  rename(`Group` = GROUP, `Age (weeks)` = AGE_WEEKS) %>%
  arrange(factor(TREATMENT, levels = sort_levels)) %>%
  select(-TREATMENT)

tab_fish <- tab_in %>%
  transmute(
    `Fish ID` = sub("grz-AD_(.*)_E", "\\1", INDIVIDUAL),
    #`Hatch date` = as.character(`Hatch date`),
    #`Sacrifice date` = as.character(`Sacrifice date`),
    #`Age at death (days)` = format(`Age (days)`, digits = 0),
    `Age at death (weeks)` = format(AGE_WEEKS, digits = 0),
    `Treatment group` = GROUP,
    #`Weight at death (g)` = DEATH_WEIGHT,
    `RNA Integrity Number (RIN)` = RIN,
    `Date of library prep` = ifelse(!is.na(LP_DATE), as.character(LP_DATE), "--"),
    `Sequenced?` = ifelse(LP > 0, "Yes", "No"),
    `Reason for exclusion` = ifelse(!is.na(REASON_FOR_EXCLUSION),
                                    REASON_FOR_EXCLUSION, "--")
  )

#------------------------------------------------------------------------------
# SAVE LATEX TABLES
#------------------------------------------------------------------------------

savetab(tab_fish, outfile_fish)
savetab(tab_summary, outfile_summary)

