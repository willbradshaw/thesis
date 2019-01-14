###############################################################################
## TABLE                                                                     ##
## Nothobranchius furzeri V/D/J co-ordinate tables                           ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")

# Configure input paths to raw tables
vh_tab_path <- "../_Data/segments/nfu/nfu_vh_tab.tsv"
jh_tab_path <- "../_Data/segments/nfu/nfu_jh_tab.tsv"
dh_tab_path <- "../_Data/segments/nfu/nfu_dh_tab.tsv"

# Configure output paths
filename_vh <- "nfu-vh-coords"
filename_jh_rss <- "nfu-jh-coords-rss"
filename_jh_seg <- "nfu-jh-coords-seg"
filename_dh_rss <- "nfu-dh-coords-rss"
filename_dh_seg <- "nfu-dh-coords-seg"

#------------------------------------------------------------------------------
# IMPORT TABLES
#------------------------------------------------------------------------------

vh_tab <- read_tsv(vh_tab_path)
jh_tab <- read_tsv(jh_tab_path)
dh_tab <- read_tsv(dh_tab_path)

#------------------------------------------------------------------------------
# SPLIT DH AND JH TABLES
#------------------------------------------------------------------------------

# Split DH into RSS info and segment info
dh_tab_seg <- dh_tab %>% select(Name, Start:End)
dh_tab_rss <- dh_tab %>% select(-(Start:End))

# Split JH into RSS info and segment info
jh_tab_seg <- jh_tab %>% select(Name, Start:End)
jh_tab_rss <- jh_tab %>% select(-(Start:End))

#------------------------------------------------------------------------------
# GENERATE LATEX TABLE
#------------------------------------------------------------------------------

savetab(vh_tab, filename_vh)
savetab(jh_tab_seg, filename_jh_seg)
savetab(jh_tab_rss, filename_jh_rss)
savetab(dh_tab_seg, filename_dh_seg)
savetab(dh_tab_rss, filename_dh_rss)
