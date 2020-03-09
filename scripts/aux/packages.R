###############################################################################
## AUX FILE                                                                  ##
## Packages                                                                  ##
###############################################################################

# List desired packages
packages <- c("plyr", # Database manipulation
              "xtable", # Export tables to LaTeX
              "Gviz", # Build genome tracks
              "GenomicRanges", # Sequence range processing
              "Biostrings", # FASTA file reading and manipulation
              "tidyverse", # Tibbles, ggplot2, dplyr etc.
              "extrafont", # Support for additional fonts
              "lubridate", # Improved date support
              "gridGraphics", # graphics-to-grid figure manipulation
              "gridExtra", # Tile grid graphics with grid.arrange
              "reshape2", # Melt tibbles and data.frames
              "BSgenome", # Extra sequence-manipulation functionality
              "DECIPHER", # Synteny analysis
              #"ggdendro", # Improved dendrogram visualisation with ggplot2
              "cowplot", # More themes and annotation functionality for ggplot2
              "ggseqlogo", # Sequence logos with ggplot2
              "ape", # Phylogenetics functionality
              "tidytree", # Tidyverse/phylo integration
              "ggtree", # Plotting phylo trees with ggplot2
              #"genoPlotR", # For chromosome synteny plot
              "DescTools", # For colour mixing
              "stringi", # String manipulations
              "alakazam", # Change-O functions 1
              #"rdi", # Repertoire Dissimilarity Index
              "shazam", # Change-O functions 2
              #"survival", # Survival curves for intro
              #"survminer" # More survival curves
              "rlang" # := helper
              )

# Define auxiliary functions
silentRequire <- function(packageName){
  # Load a package silently if possible
  suppressMessages(
    require(packageName, character.only = TRUE, quietly = TRUE)
  )
}

# Load/install packages as appropriate
for (p in packages) silentRequire(p)

# Install RDI from local source if needed
if ("rdi" %in% names(snakemake@params)){
    suppressMessages(install.packages(snakemake@params[["rdi"]],
                                      repos=NULL, type="source"))
    silentRequire("rdi")
}
