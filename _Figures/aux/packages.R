###############################################################################
## AUX FILE                                                                  ##
## Packages                                                                  ##
###############################################################################

# Make sure installation packages are available
if (!("BiocManager" %in% installed.packages())) install.packages("BiocManager")
if (!("devtools" %in% installed.packages())) install.packages("devtools")

# List desired packages
packages <- c("Gviz", # Build genome tracks
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
              "ggdendro", # Improved dendrogram visualisation with ggplot2
              "cowplot", # More themes and annotation functionality for ggplot2
              "ggseqlogo", # Sequence logos with ggplot2
              "ape", # Phylogenetics functionality
              "tidytree", # Tidyverse/phylo integration
              "ggtree" # Plotting phylo trees with ggplot2
              )

# Define auxiliary functions
silentRequire <- function(packageName){
  # Load a package silently if possible
  suppressMessages(
    require(packageName, character.only = TRUE, quietly = TRUE)
  )
}
quietPackage <- function(packageName){
  # Load a package silently if possible, otherwise install it quietly
  if (!silentRequire(packageName)){
    BiocManager::install(packageName, dep = TRUE, quietly = TRUE)
    if(!silentRequire(packageName)){
      stop("Could not install package", packageName) 
    }
  }
}

# Load/install packages as appropriate
for (p in packages) quietPackage(p)

# Clear namespace
rm(list = c("p", "packages"))