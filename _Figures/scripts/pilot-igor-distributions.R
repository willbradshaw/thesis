###############################################################################
## FIGURE                                                                    ##
## IGOR PLOTS FOR PILOT DATASET                                              ##
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

# Specify input paths
experiment <- "pilot"
segments_path_individual <- paste0("../_Data/igor/segments/", experiment,
                                   "-individual-segments.tsv")
segments_path_group <- paste0("../_Data/igor/segments/", experiment,
                                   "-group-segments.tsv")
indels_path_individual <- paste0("../_Data/igor/indels/", experiment,
                                   "-individual-indels.tsv")
indels_path_group <- paste0("../_Data/igor/indels/", experiment,
                              "-group-indels.tsv")
v_names_path <- "../_Data/segments/nfu/nfu_vh_name_conversion.csv"
d_names_path <- "../_Data/segments/nfu/nfu_dh_name_conversion.csv"
j_names_path <- "../_Data/segments/nfu/nfu_jh_name_conversion.csv"

# Specify parameters
individuals <- paste0("2-0", seq(3, 6))
palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
             colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]])
# Output path
filename_base <- paste0(experiment, "-igor")

#------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
#------------------------------------------------------------------------------

import_segments <- function(path, v_names, d_names, j_names){
  col <- cols(p = "d", .default = "c")
  tab <- import_tsv(path, col = col) %>%
    mutate(SEGMENT_TYPES = paste0(ifelse(is.na(v), "", "V"),
                                  ifelse(is.na(d), "", "D"),
                                  ifelse(is.na(j), "", "J"))) %>%
    rename(ID = id, P = p, V_NAME_NEW = v, D_NAME_NEW = d, J_NAME_NEW = j) %>%
    full_join(v_names, by = "V_NAME_NEW") %>%
    full_join(d_names, by = "D_NAME_NEW") %>%
    full_join(j_names, by = "J_NAME_NEW") %>%
    mutate(V_NAME_OLD = gsub("_", " / ", V_NAME_OLD),
           D_NAME_OLD = gsub("_", " / ", D_NAME_OLD),
           J_NAME_OLD = gsub("_", " / ", J_NAME_OLD))
}

plot_segments <- function(segments, grouped = FALSE){
  indiv <- segments_individual %>% filter(SEGMENT_TYPES == segments)
  g_base <- ggplot(mapping = aes_string(x=paste0(segments, "_NAME_OLD"), y="P",
                                        colour="ID", group="ID")) +
    scale_colour_manual(values = palette, name = "Individual") +
    scale_y_continuous(name = "Probability (%)",
                       labels = function(y) y * 100,
                       expand=c(0,0)) +
    theme_classic() + theme_base + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = fontsize_base * fontscale_title,
                                 angle = 45, hjust = 1),
      legend.title = element_text(margin = margin(r = 0.5, unit = "cm"))
    )
  if (!grouped) return(g_base + geom_line(data = indiv))
  group <- segments_group %>% filter(SEGMENT_TYPES == segments)
  g <- g_base + geom_line(data = indiv, alpha = 0.8, linetype = 2) +
    geom_line(data = group, size = 1.5, 
              colour = colours_igseq[["ageing_group2"]])
  return(g)
}

plot_indels <- function(ev, grouped = FALSE){
  indiv <- indels_individual %>% filter(event == ev)
  g_base <- ggplot(mapping = aes(x=n, y=p, colour=id, group=id)) +
    scale_colour_manual(values = palette, name = "Individual") +
    scale_y_continuous(name = "Probability (%)",
                       labels = function(y) y * 100, expand = c(0,0)) +
    theme_classic() + theme_base + theme(
      legend.title = element_text(margin = margin(r = 0.5, unit = "cm"))
    )
  if (!grouped) return(g_base + geom_line(data = indiv))
  group <- indels_group %>% filter(event == ev)
  g <- g_base + geom_line(data = indiv, alpha = 0.8, linetype = 2) +
    geom_line(data = group, size = 1.5, 
              colour = colours_igseq[["ageing_group2"]])
  return(g)
}

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

# Segment name correspondence
v_names <- suppressMessages(read_csv(v_names_path))
d_names <- suppressMessages(read_csv(d_names_path))
j_names <- suppressMessages(read_csv(j_names_path))

# Segments
segments_individual <- import_segments(segments_path_individual, v_names,
                                       d_names, j_names)
segments_group <- import_segments(segments_path_group, v_names, d_names,
                                  j_names)

# Indels
col_indels <- cols(p = "d", n = "i", .default = "c")
indels_individual <- import_tsv(indels_path_individual, col_indels)
indels_group <- import_tsv(indels_path_group, col_indels)

#------------------------------------------------------------------------------
# VISUALISE (AND SAVE) SEGMENT AND INDEL DISTRIBUTIONS
#------------------------------------------------------------------------------
C_Appendices/Tables
# Segment usage plot
g_d <- plot_segments("D", TRUE)
g_j <- plot_segments("J", TRUE)
g_dj <- gplot_grid(g_d + theme(legend.position = "none"), 
                   g_j + theme(legend.position = "none"), 
                   ncol = 2, nrow = 1, labels = c("B", "C"))
g_v <- plot_segments("V", TRUE)
g_segments <- gplot_grid_onelegend(g_v, g_dj, nrow = 2, ncol = 1, 
                                   labels = c("A", ""), plot_width = 30,
                                   plot_height = 25)

# Indel plot
xcont <- function(title, max) scale_x_continuous(name = title, expand = c(0,0),
                                                 limits = c(NA,max))
g_ins_vd <- plot_indels("vd_ins", TRUE) + xcont("# V/D-insertions", 20)
g_ins_dj <- plot_indels("dj_ins", TRUE) + xcont("# V/D-insertions", 20)
g_del_v <- plot_indels("v_3_del", TRUE) + xcont("# 3' V-deletions", 20)
g_del_j <- plot_indels("j_5_del", TRUE) + xcont("# 5' J-deletions", 20)
g_del_d3 <- plot_indels("d_3_del", TRUE) + xcont("# 3' D-deletions", 25)
g_del_d5 <- plot_indels("d_5_del", TRUE) + xcont("# 5' D-deletions", 25)
g_indels <- gplot_grid_onelegend(g_ins_vd, g_ins_dj, g_del_v, g_del_d5, 
                                 g_del_d3, g_del_j, nrow = 3, ncol = 2, 
                                 plot_width = 21,  plot_height = 29.7)

savefig(g_segments, paste0(filename_base, "-segments"),
        height = 25, width = 30)
savefig(g_indels, paste0(filename_base, "-indels"),
        height = 29.7, width = 21)