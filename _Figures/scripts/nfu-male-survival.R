###############################################################################
## FIGURE                                                                    ##
## N. furzeri GRZ survival curves                                            ##
###############################################################################

#------------------------------------------------------------------------------
# PREAMBLE
#------------------------------------------------------------------------------

# Source auxiliary files (packages, fonts, etc.)
source("aux/packages.R")
source("aux/fonts.R")
source("aux/palette.R")
source("aux/io.R")
source("aux/ggplot2.R")

# Configure input paths
survival_path <- "../_Data/fish_data/grz-survival.txt"

# Configure parameters
group_size <- 28
start_date = "2016-12-12"

# Output options
plot_height <- 18
plot_ratio <- 0.9
filename <- "nfu-male-survival"

#------------------------------------------------------------------------------
# IMPORT AND CONFIGURE DATA
#------------------------------------------------------------------------------

# Import data
survival <- suppressMessages(read_tsv(survival_path)) %>%
  filter(type < 3) %>% # Filter to male fish
  mutate(time_weeks = time / 7) # Convert days to weeks

#------------------------------------------------------------------------------
# COMPUTE SURVIVAL CURVE
#------------------------------------------------------------------------------

# Make survival object
S <- Surv(survival[["time"]], survival[["event"]])

# Fit survival curve (no subgroups)
SF <- survfit(S ~ 1, data = survival)

#------------------------------------------------------------------------------
# PLOT SURVIVAL
#------------------------------------------------------------------------------

# Get location of median lines
x_median <- SF$time[min(which(SF$surv < 0.5))]
y_median <- 0.5
abline_data <- data.frame(x=x_median, y=y_median)

# Survplot objects
col <- "#66A61E"
g_surv <- ggsurvplot(SF, size = 2, conf.int = FALSE, palette = col)$plot +
  geom_segment(data=abline_data, aes(x=0, xend=x, y=y, yend=y), lty=2, size=1, colour=col) +
  geom_segment(data=abline_data, aes(x=x, xend=x, y=0, yend=y), lty=2, size=1, colour=col) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1.05), labels = function(y) y * 100,
                     name = "Survival probability (%)") +
  scale_x_continuous(expand=c(0,0), limits=c(0,NA), name = "Age (weeks)",
                     labels = function(x) x / 7, breaks = function(a) seq(a[1],a[2],35),
                     sec.axis = sec_axis(trans = ~., name = "Age (days)")) +
  scale_colour_manual(name = "", labels = "Male GRZ-AD\n(n = 52)", values = col) +
  theme_classic() + theme_base + theme(
    axis.title.x.top = element_text(margin = margin(b = 0.6, t = 0, unit = "cm"),
                                    size = fontsize_base * fontscale_label),
    axis.title.x.bottom = element_text(size = fontsize_base * fontscale_label),
    legend.position = c(0.7, 0.85),
    plot.title = element_blank(),
    plot.margin = margin(t = 0.4, b = 0.4, r = 0.8, l = 0.4, unit = "cm"),
    legend.text = element_text(size = fontsize_base * fontscale_title),
    axis.text = element_text(size = fontsize_base * fontscale_title),
    axis.text.x.bottom = element_text(margin = margin(t=0.2,unit = "cm")),
    axis.text.x.top = element_text(margin = margin(b=0.2,unit = "cm")),
    axis.text.y = element_text(margin = margin(r = 0.1, unit = "cm")),
    axis.title.y = element_text(size = fontsize_base * fontscale_label)
  )


#------------------------------------------------------------------------------
# SAVE PLOT
#------------------------------------------------------------------------------

savefig(plot = g_surv, filename = filename,
        height = plot_height, ratio = plot_ratio)
