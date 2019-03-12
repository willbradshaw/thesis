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
entropies_path_individual <- paste0("../_Data/igor/entropies/", experiment,
                                   "-individual-entropies.tsv")
entropies_path_group <- paste0("../_Data/igor/entropies/", experiment,
                              "-group-entropies.tsv")
entropies_path_individual <- "../_Data/igor/entropies/pilot-individual-entropies.tsv"
entropies_path_group <- "../_Data/igor/entropies/pilot-group-entropies.tsv"
v_names_path <- "../_Data/segments/nfu/nfu_vh_name_conversion.csv"
d_names_path <- "../_Data/segments/nfu/nfu_dh_name_conversion.csv"
j_names_path <- "../_Data/segments/nfu/nfu_jh_name_conversion.csv"

# Specify parameters
individuals <- paste0("2-0", seq(3, 6))
palette <- c(colours_igseq[["pilot_rep1"]], colours_igseq[["pilot_rep2"]],
             colours_igseq[["pilot_rep3"]], colours_igseq[["pilot_rep4"]])
entropy_events <- tibble(event = c("model", "v_choice", "d_gene", "j_choice", 
                                  "vd_dinucl", "vd_ins", "dj_dinucl", "dj_ins",
                                  "v_3_del", "d_5_del", "d_3_del", "j_5_del"),
                         category = c("Recombination events (total)", rep("Gene choice", 3), 
                                      rep("Insertions",4),rep("Deletions",4)),
                         short = c("total", "V", "D", "J", 
                                   "VD nts", "VD\nlen", "DJ nts", "DJ\nlen",
                                   "delV", "delD", "delD", "delJ")
)

# Output path
filename_base <- "pilot-igor"

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
                       labels = function(y) y * 100) +
    theme_classic() + theme_base + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = fontsize_base * fontscale_title,
                                 angle = 45, hjust = 1)
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
                       labels = function(y) y * 100) +
    theme_classic() + theme_base
  if (!grouped) return(g_base + geom_line(data = indiv))
  group <- indels_group %>% filter(event == ev)
  g <- g_base + geom_line(data = indiv, alpha = 0.8, linetype = 2) +
    geom_line(data = group, size = 1.5, 
              colour = colours_igseq[["ageing_group2"]])
  return(g)
}

import_entropies <- function(path, ...){
  col_entropies <- cols(h = "d", .default = "c")
  htab <- import_tsv(path, col_entropies) %>% 
    full_join(entropy_events, by = "event") %>%
    group_by(category, short) %>% summarise(h = sum(h)) # ...
  return(htab)
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

# Entropies
entropies_individual <- import_entropies(entropies_path_individual)
entropies_group <- import_entropies(entropies_path_group)

#------------------------------------------------------------------------------
# VISUALISE (AND SAVE) SEGMENT AND INDEL DISTRIBUTIONS
#------------------------------------------------------------------------------

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
g_ins_vd <- plot_indels("vd_ins", TRUE) + xlab("# V/D-insertions") + xlim(c(NA,20))
g_ins_dj <- plot_indels("dj_ins", TRUE) + xlab("# D/J-insertions") + xlim(c(NA,20))
g_del_v <- plot_indels("v_3_del", TRUE) + xlab("# 3' V-deletions")
g_del_j <- plot_indels("j_5_del", TRUE) + xlab("# 5' J-deletions")
g_del_d3 <- plot_indels("d_3_del", TRUE) + xlab("# 3' D-deletions") + xlim(c(NA,25))
g_del_d5 <- plot_indels("d_5_del", TRUE) + xlab("# 5' D-deletions") + xlim(c(NA,25))
g_indels <- gplot_grid_onelegend(g_ins_vd, g_ins_dj, g_del_v, g_del_d5, 
                                 g_del_d3, g_del_j, nrow = 3, ncol = 2, 
                                 plot_width = 21,  plot_height = 29.7)

savefig(g_segments, paste0(filename_base, "-segments"),
        height = 25, width = 30)
savefig(g_indels, paste0(filename_base, "-indels"),
        height = 29.7, width = 21)

# TODO: Remove gap between axis and 0 height/width

#------------------------------------------------------------------------------
# VISUALISE (AND SAVE) SEGMENT AND INDEL DISTRIBUTIONS
#------------------------------------------------------------------------------

# Obtain entropy composition of grouped model
categories <- entropy_events %>% pull(category) %>% unique
shorts <- entropy_events %>% pull(short) %>% unique %>% 
  (function(x) c(categories, x[-1]))
# shorts <- c(categories, "V", "D", "J", "VD nts", "VD\nlength", "DJ nts", 
#             "DJ\nlength", "delV", "delD", "delJ")
alphas <- c(1,0.7,0.7,0.7,1,0.4,1,1,0.4,1,0.4,1,0.4,1)
h_total_group <- entropies_group %>% filter(short == "total") %>% pull(h)

h_group <- entropies_group %>% group_by(category) %>% 
  summarise(h = sum(h)) %>% mutate(short = category) %>% 
  bind_rows(entropies_group %>% filter(short != "total")) %>%
  mutate(level = ifelse(category == "Recombination events (total)", 1, 
                        ifelse(short == category, 2, 3))) %>%
  mutate(label = ifelse(level == 3, short,
                        paste0(short, ": ", round(h, 2), " bits")),
         short = factor(short, levels = shorts),
         category = factor(category, levels = categories)) %>%
  arrange(short) %>% group_by(level) %>%
  mutate(ch = cumsum(h), chl = lag(ch, default = 0),
         chx = (ch + chl)/2, chx = h_total_group - chx)

# Plot entropy composition
g_h_group <- ggplot(h_group) +
  geom_col(aes(x=level, y=h, fill=category, alpha=short), width = 1) +
  geom_text(aes(x=level, y = chx, label = label),
            hjust = 0.5, vjust = 0.5, family = font, size = 5.5) +
  scale_alpha_manual(values = alphas) +
  coord_flip() + scale_x_reverse() + scale_y_reverse() +
  theme_minimal() + theme_base + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.position = "none"
  )
