###############################################################################
## FIGURE                                                                    ##
## IGOR ENTROPY COMPOSITION FOR AGEING DATASET                               ##
###############################################################################

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "aux.R"))
parse_snakemake(snakemake)

write_log("Parsed global Snakemake properties.")
write_log("Loaded packages and auxiliary functions.")

# Specify parameters
age_groups <- c("39", "56", "73", "128")
entropy_events <- tibble(event = c("model", "v_choice", "d_gene", "j_choice", 
                                  "vd_dinucl", "vd_ins", "dj_dinucl", "dj_ins",
                                  "v_3_del", "d_5_del", "d_3_del", "j_5_del"),
                         category = c("Recombination events (total)", rep("Gene choice", 3), 
                                      rep("Insertions",4),rep("Deletions",4)),
                         short = c("total", "V", "D", "J", 
                                   "VD nts", "VD\nlen", "DJ nts", "DJ\nlen",
                                   "delV", "delD", "delD", "delJ")
)
palette <- c("#FFCCCC","#66B266","#8080E6","#E66666")

#------------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
#------------------------------------------------------------------------------

import_entropies <- function(path){
  col_entropies <- cols(h = "d", .default = "c")
  htab <- import_tsv(path, col_entropies) %>% 
    full_join(entropy_events, by = "event") %>%
    group_by(id, category, short) %>% summarise(h = sum(h))
  return(htab)
}

entropy_table <- function(entropies, events = entropy_events){
  # Get levels and references from events table
  categories <- events %>% pull(category) %>% unique
  shorts <- events %>% pull(short) %>% unique %>% 
    (function(x) c(categories, x[-1]))
  top_short <- events %>% filter(event == "model") %>% pull(short)
  top_cat <- events %>% filter(event == "model") %>% pull(category)
  h_total <- entropies %>% filter(short == "total") %>% pull(h)
  # Construct table
  htab <- entropies_group %>% group_by(id, category) %>% 
    summarise(h = sum(h)) %>% mutate(short = category) %>% 
    bind_rows(entropies_group %>% filter(short != top_short)) %>%
    mutate(level = ifelse(category == top_cat, 1, 
                          ifelse(short == category, 2, 3))) %>%
    mutate(label = ifelse(level == 3, short,
                          paste0(short, ": ", round(h, 2), " bits")),
           short = factor(short, levels = shorts),
           category = factor(category, levels = categories)) %>%
    arrange(short) %>% group_by(id, level) %>%
    mutate(ch = cumsum(h), chl = lag(ch, default = 0),
           chx = (ch + chl)/2, chx = sum(h) - chx) %>%
    group_by(id, category) %>% select(-ch, -chl) %>%
    mutate(alpha = ifelse(level == 1, 1, ifelse(level == 2, 0.7, c(0.4,1))))
  return(htab)
}

plot_entropy_table <- function(htab, pal = palette){
  # Create stacked entropy composition plot from prepared entropy table
  g <- ggplot(htab) + theme_minimal() + theme_base + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    plot.margin = margin(t = 1.2, l = -1, r = -1, unit = "cm"),
    #panel.border = element_blank(), panel.background = element_blank(), 
    legend.position = "none"
  )
  # Create basic plot structure
  g <- g + geom_col(aes(x=level, y=h, fill=category, alpha=short), width = 1) +
    geom_text(aes(x=level, y = chx, label = label),
              hjust = 0.5, vjust = 0.5, family = font, size = 5.5) +
    coord_flip() + scale_x_reverse() + scale_y_reverse()
  
  # Configure colour display
  g <- g + scale_alpha_manual(values = h_group %>% arrange(id, level) %>% 
                                pull(alpha)) +
    scale_fill_manual(values = pal)
  return(g)
}

#------------------------------------------------------------------------------
# IMPORT DATA
#------------------------------------------------------------------------------

entropies_group <- import_entropies(inpath_entropies_group) %>%
  mutate(age_days = factor(id, levels = age_groups))

entropies_individual <- import_entropies(inpath_entropies_indiv) %>%
  mutate(age_group = as.integer(sub("-.*", "", id)),
         age_days = age_groups[age_group])

#------------------------------------------------------------------------------
# VISUALISE ENTROPY COMPOSITION
#------------------------------------------------------------------------------

# Process entropies for plotting
h_group <- entropy_table(entropies_group) %>% 
  mutate(age_days = factor(id, levels = age_groups)) %>%
  filter(id != "128")

# Plot entropy composition
g_h_group <- plot_entropy_table(h_group) + 
  facet_wrap(~id, ncol = 1, 
             labeller = as_labeller(function(c) paste0("Age (days) = ", c))) +
  theme(strip.text.x = element_text(hjust = 0.9))
ggsave(plot = g_h_group, filename = outpath_plot, device = "svg", units = "cm",
       height=20, width = 28)
