###############################################################################
## AUX FILE                                                                  ##
## Auxiliary functions for tree processing                                   ##
###############################################################################

source(file.path(aux_dir, "ggplot2.R"))

# Check tip labelling and find ancestry patterns
grepl_tips <- function(pattern, tree_object){
  # Search for a string pattern in tip labels and return all matching tips
  tree_object$tip.label[grepl(pattern, tree_object$tip.label)]
}

get_mrca <- function(pattern, tree_object){
  # Find the MRCA of all tips matching a string pattern
  nodes <- grepl_tips(pattern, tree_object)
  if (length(nodes) <= 1) return(NA)
  return(MRCA(tree_object, nodes))
}

check_monophyly <- function(pattern, tree_object, 
                            internal_numeric = TRUE){
  # Check that nodes matching a tip pattern form a monophyletic group
  mrca <- get_mrca(pattern, tree_object)
  clade <- offspring(.data = tree_object, .node = mrca)
  subtree <- as_tibble(tree_object) %>% filter(node %in% clade) %>%
    mutate(label = ifelse(label == "", 0, label))
  if (internal_numeric){ # Internal nodes marked with numeric support values; tip labels non-numeric
    subtree <- suppressWarnings(subtree %>% filter(is.na(as.integer(label))))
  } else {
    subtree <- subtree %>% filter(!is.na(label))
  }
  return(ifelse(all(grepl(pattern, subtree$label)), mrca, FALSE))
}

#------------------------------------------------------------------------------
# Process tree tables
#------------------------------------------------------------------------------

label_tips <- function(tree_tab, internal_numeric = TRUE){
  # Identify and label tip nodes in a tree tibble
  if (internal_numeric){
    tree_tab[["tip"]] <- suppressWarnings(is.na(as.numeric(tree_tab[["label"]])))
  } else {
    tree_tab[["tip"]] <- !is.na(tree_tab[["label"]])
  }
  tree_tab[["tip"]][which(tree_tab[["label"]] == "")] <- FALSE
  return(tree_tab)
}

column_from_label <- function(tree_tab, column, pattern, replacement,
                              internal_numeric = TRUE){
  # Infer a new column from tip labels, leaving internal nodes blank
  if (!"tip" %in% colnames(tree_tab)){
    tab <- label_tips(tree_tab)
  } else {
    tab <- tree_tab
  }
  tab[[column]] <- sub(pattern, replacement, tab[["label"]])
  tab[[column]] <- ifelse(tab[["tip"]], tab[[column]], NA)
  return(tab)
}

propagate_column <- function(tree_tab, column){
  # Propagate column values to internal nodes
  nodes <- tree_tab[is.na(tree_tab[[column]]),] %>% pull(node)
  node_depth <- sapply(nodes, function(n) offspring(tree_tab, n) %>% nrow)
  nodes <- nodes[order(node_depth)] # Sort by fewest total offspring
  while (length(nodes) > 0){
    for (n in nodes){
      row <- which(tree_tab[["node"]] == n)
      if (is.na(tree_tab[[column]][row])){
        subtree <- offspring(tree_tab, n)
        subcol <- unique(subtree[[column]])
        if (length(subcol) == 1){
          tree_tab[[column]][row] <- subcol[1]
        } else if ("CONFLICT" %in% subcol) {
          tree_tab[[column]][row] <- "CONFLICT"
        } else if (sum(!is.na(subcol)) > 1){
          tree_tab[[column]][row] <- "CONFLICT"
        } else {
          tree_tab[[column]][row] <- NA
        }
      }
    }
    nodes <- tree_tab[is.na(tree_tab[[column]]),] %>% pull(node)
    node_depth <- sapply(nodes, function(n) offspring(tree_tab, n) %>% nrow)
    nodes <- nodes[order(node_depth)] # Sort by fewest total offspring
  }
  # TODO: Convert CONFLICT nodes back to NAs?
  tree_tab[tree_tab[[column]] == "CONFLICT",][[column]] <- NA
  return(tree_tab)
}

#------------------------------------------------------------------------------
# Default tree theme
#------------------------------------------------------------------------------

theme_tre <- theme_tree2() + theme_base +
  theme(legend.position = "top",
        legend.margin = margin(b=0,t=0.5, unit="cm"),
        plot.margin = margin(t=0, unit="cm"),
        axis.title.x = element_text(margin = margin(t=0, b=0.5, unit="cm")),
        axis.title.y = element_text(margin = margin(l = 2, t = 2, unit = "cm")))
