###############################################################################
## AUX FILE                                                                  ##
## Font settings                                                             ##
###############################################################################

# Load extra fonts
suppressMessages(font_install("fontcm", prompt=FALSE))
suppressMessages(loadfonts(device="postscript"))
suppressMessages(loadfonts(device="pdf"))

# Configure font settings
font <- "CM Sans" # Main figure font
titlefont <- font # Font for axis titles etc. (if different)
fontsize_base <- 12 # Basic figure font size
fontscale_title <- 1.5 # Default axis-title scale relative to regular font
fontscale_main <- 1.5 # Default plot-title scale
fontscale_label <- 2.5 # Default subfigure-label scale (A, B, etc)
