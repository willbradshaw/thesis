###############################################################################
## TABLE                                                                     ##
## Ageing cohort fish and group summary                                      ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/io.R")

# Path to input data
inpath <- "../_Data/fish_data/igseq-cohorts.csv"

# Output paths
outfile_fish <- "igseq-cohorts-fish"
outfile_summary <- "igseq-cohorts-summary"

#------------------------------------------------------------------------------
# IMPORT DATA AND GENERATE SUMMARY
#------------------------------------------------------------------------------

tab_fish <- suppressMessages(read_csv(inpath))

tab_summary <- tab_fish %>% group_by(Group) %>%
  summarise(`# Fish` = n(), `Hatch date` = as.character(first(`Hatch date`)),
            `Sacrifice date` = as.character(first(`Sacrifice date`)),
            `Age (days)` = format(first(`Age (days)`), digits = 0),
            `Age (weeks)` = format(first(`Age (weeks)`), digits = 3),
            `Mean weight (g)` = format(mean(`Death weight (g)`), digits = 3)
            )

tab_fish <- tab_fish %>%
  mutate(`Hatch date` = as.character(`Hatch date`),
         `Sacrifice date` = as.character(`Sacrifice date`),
         `Age (days)` = format(`Age (days)`, digits = 0),
         `Age (weeks)` = format(`Age (weeks)`, digits = 3),
         `Fish ID` = sub("grz-AD_(.*)_E", "\\1", `Fish ID`)
  ) %>% select(-Strain, -Sex)


#------------------------------------------------------------------------------
# SAVE LATEX TABLES
#------------------------------------------------------------------------------

savetab(tab_fish, outfile_fish)
savetab(tab_summary, outfile_summary)

