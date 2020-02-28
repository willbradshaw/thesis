# AUXILIARY R FUNCTIONS FOR SNAKEMAKE SCRIPTS

aux_dir <- snakemake@params[["aux"]]
source(file.path(aux_dir, "packages.R"))
source(file.path(aux_dir, "snakemake.R"))
source(file.path(aux_dir, "blast.R"))
source(file.path(aux_dir, "fonts.R"))
source(file.path(aux_dir, "palette.R"))
source(file.path(aux_dir, "ggplot2.R"))
source(file.path(aux_dir, "gviz.R"))
source(file.path(aux_dir, "io.R"))
source(file.path(aux_dir, "ranges.R"))
source(file.path(aux_dir, "trees.R"))
source(file.path(aux_dir, "changeo.R"))
