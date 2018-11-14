###############################################################################
## AUX FILE                                                                  ##
## Configure I/O defaults                                                    ##
###############################################################################

# Source packages
source("aux/packages.R")

# Prepare grid layout for multi-figure plotting
split_layout <- function(width, height, width_ratios = 1, 
                         height_ratios = 1, unit = "cm"){
  fig_heights <- sapply(height_ratios, 
                        function(x) x/sum(height_ratios) * height)
  fig_widths <- sapply(width_ratios, 
                       function(x) x/sum(width_ratios) * width)
  layout <- grid.layout(
    ncol = length(fig_widths),
    nrow = length(fig_heights),
    heights = unit(fig_heights, plot_unit),
    widths = unit(fig_widths, plot_unit)
  )
  return(layout)
}

# ggsave function for grid objects
savefig <- function(plot, filename, device, height = NA, 
                    width = NA, ratio = NA,
                    units = "cm", dpi = 320){
  if (!dir.exists(device)) dir.create(device, recursive = TRUE)
  ggsave(plot = plot,
         filename = file.path(device, paste0(filename, ".", device)),
         device = device,
         height = height,
         width = ifelse(is.na(width), height * ratio, width),
         units = units,
         dpi = dpi
  )
}

