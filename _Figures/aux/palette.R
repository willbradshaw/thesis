###############################################################################
## AUX FILE                                                                  ##
## Colour palette settings                                                   ##
###############################################################################

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")

# Specify font colours
fontcolour_dark <- "black"
fontcolour_mid <- "#606060"
fontcolour_alert <- "red"

# Specify colour palette for locus chapter
colours <- list(
  VH = "#E41A1C",
  DH = "#377EB8",
  JH = "#4DAF4A",
  CD = "#984EA3",
  CM = "#FF7F00",
  X1 = "#FFFF33",
  X2 = "#A65628",
  CZ = "#F781BF",
  GR = "#444444", # Grey
  HL1 = "#e31a1c",
  HL2 = "#1f78b4",
  WH = "#FFFFFF", # White
  LG = "#DDDDDD", # Light grey
  LR = "#FFCCCC", # Light red
  RD = "#E31919", # Red
  nfu = "#66c2a5",
  xma = "#fc8d62",
  ola = "#8da0cb"
)


# Specify colour palette for igseq chapter
colours_igseq <- list(
  exp_pilot = "#1F78B4",
  exp_ageing = "#E31A1C",
  exp_gut = "#33A02C",
  pilot_rep1 = "#999999",
  pilot_rep2 = "#E69F00",
  pilot_rep3 = "#56B4E9",
  pilot_rep4 = "#009E73",
  ageing_group1 = "#1b9e77", # 39 days
  ageing_group2 = "#d95f02", # 56 days
  ageing_group3 = "#7570b3", # 73 days
  ageing_group4 = "#e7298a", # 128 days
  gut_group1 = "#e78ac3", # YI_6wk
  gut_group2 = "#8da0cb", # WT_16wk
  gut_group3 = "#e5c494", # ABX_16wk
  gut_group4 = "#fc8d62", # SMT_16wk
  gut_group5 = "#a6d854", # YMT_16wk
  #gut_allold = "#b3b3b3" # All 16-week groups
  gut_allold = "#ffd92f" # All 16-week groups
  
)