###############################################################################
## AUX FILE                                                                  ##
## Configure ggplot2 display defaults                                        ##
###############################################################################

# Auxiliary palette function
gg_color_hue <- function(n) {
  # Get default ggplot2 colour palette for n colours as a character vector
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Define base theme
theme_base <-   theme(
  legend.position = "bottom",
  axis.text = element_text(size = fontsize_base, family = font, 
                           colour = "black"),
  axis.title.y = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(r=3,unit="mm")),
  axis.title.x = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t=5,unit="mm")),
  legend.text = element_text(size = fontsize_base, 
                             family = font, colour = "black"),
  legend.title = element_text(size = fontsize_base,
                              family = font, colour = "black",
                              face = "bold", vjust=0.5, 
                              margin=margin(r=3, unit="mm")),
  plot.title = element_text(size = fontsize_base * fontscale_main,
                            family = titlefont, colour="black",
                            hjust = 0.5, vjust = 0.5,
                            margin = margin(b=5, unit="mm")),
  plot.margin = margin(t=1.5, l=0.5, r=0.5, b = 0.5, unit="cm"),
  strip.background = element_blank(),
  strip.text.x = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(t = 3, b = 3, unit="mm")),
  strip.text.y = element_text(size = fontsize_base * fontscale_title,
                              family = titlefont, colour = "black",
                              margin = margin(l=3, r=3, unit="mm")),
  legend.justification = "center"
)

# Combine plots into a grid with plot_grid
gplot_grid <- function(..., ncol = 2, nrow = 1, labels = "AUTO"){
  plot_grid(..., ncol = ncol, nrow = nrow, labels = labels,
            label_fontfamily = titlefont, label_fontface = "plain",
            label_size = fontsize_base*fontscale_label)
}

# Plot two plots side-by-side with a common legend
gplot_grid_onelegend <- function(..., plot_height, plot_width,
                                 ncol = 2, nrow = 1, plot_unit = "cm",
                                 labels = "AUTO"){
  plotlist <- list(...)
  # Extract legend
  g <- ggplotGrob(plotlist[[1]] + theme(legend.position = "bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  # Combine plots without legend
  plotlist <- lapply(plotlist, function(p) p + theme(legend.position = "none"))
  plt <- gplot_grid(plotlist = plotlist,
                    ncol = ncol, nrow = nrow, labels = labels)
  combined <- arrangeGrob(plt, legend, ncol = 1, nrow = 2,
                          heights = unit.c(unit(1, "npc") - lheight, lheight))
  
  # Visualise plot
  map_layout <- grid.layout(ncol = 1, nrow = 1, 
                            heights = unit(plot_height, plot_unit), 
                            widths = unit(plot_width, plot_unit)
  )
  vtop <- viewport(layout = map_layout)
  grid.newpage()
  pushViewport(vtop)
  pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
  grid.draw(combined)
  popViewport(1)
  plt_out <- grid.grab()
  return(plt_out)
}
