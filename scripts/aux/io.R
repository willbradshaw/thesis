###############################################################################
## AUX FILE                                                                  ##
## Configure I/O defaults                                                    ##
###############################################################################

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
savefig <- function(plot, filename, device = c("svg", "png"), height = NA, 
                    width = NA, ratio = NA,
                    units = "cm", dpi = 320){
  for (d in device){
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    ggsave(plot = plot,
           filename = file.path(d, paste0(filename, ".", d)),
           device = d,
           height = height,
           width = ifelse(is.na(width), height * ratio, width),
           units = units,
           dpi = dpi, limitsize = FALSE
    )
  }
}

# Save tabular object to LaTeX with xtable
savetab <- function(tab, filename, auto = TRUE, align=NULL, digits = NULL, 
                    display = NULL, hline.after = NULL, 
                    include.rownames = FALSE, include.colnames = TRUE,
                    rotate.colnames = FALSE){
  # Make xtable object
  if (auto){
    xt <- xtable(tab, auto = auto)
  } else {
    xt <- xtable(tab, align = align, digits = digits, display = display)
  }

    # Print to file
  xtable::print.xtable(xt, 
                       file = ifelse(is.null(filename), "", 
                                     file.path("tables", 
                                               paste0(filename, ".tex"))),
                       floating = FALSE,
                       hline.after = hline.after,
                       include.rownames = include.rownames,
                       include.colnames = include.colnames,
                       rotate.colnames = rotate.colnames,
                       add.to.row = list(pos = list(-1, 0, nrow(xt)),
                                         command = c('\\toprule ', 
                                                     '\\midrule ',
                                                     '\\bottomrule ')))
}

# Save a quantity as raw text
savetxt <- function(x, path, ncolumns = 1, append = FALSE, sep = ""){
  write(x, file = path, ncolumns = ncolumns, append = append, sep=sep)
}

# Import and export Change-O DBs
import_tsv <- function(path, col = cols(.default = col_character())){
  suppressMessages(read_tsv(path, col_types = col))
}
