################################################################
## AUXILIARY R FUNCTIONS FOR MANIPULATING SNAKEMAKE VARIABLES ##
################################################################

#-----------------------------------------------------------------------------
# Parsing Snakemake object
#-----------------------------------------------------------------------------

abspath <- function(path, dir=getwd()){
    # Add a directory to a relative file path to make it absolute
    return(file.path(dir, path))
}

assign_global <- function(name, value, abs=FALSE){
    # Assign a value to a name in the global environment
    # HANDLE WITH CARE! Not for normal functions - just snakemake assignment
    if (abs) value <- abspath(value)
    assign(name, value, envir = .GlobalEnv)
}

parse_named <- function(object, prefix, abs=FALSE){
    # Parse elements of an object with named elements into separate objects
    for (n in names(object)){
        name <- ifelse(nchar(prefix) > 0, paste0(prefix, "_", n), n)
        assign_global(name, object[[n]], abs)
    }
}

parse_unnamed <- function(object, prefix, abs=FALSE){
    # Parse elements of an object without named elements into separate objects
    if (length(object) == 1){
        assign_global(prefix, object[[1]], abs)
    } else {
        assign_global(paste0(prefix, "s"), object, abs)
    }
}

parse_snakemake_slot <- function(snakemake, slot_name, prefix, abs=FALSE){
    # Parse elements of a snakemake object slot into separate objects
    o <- slot(snakemake, slot_name)
    if (length(o) == 0) return()
    if (is.null(names(o))) return(parse_unnamed(o, prefix, abs))
    o <- o[str_count(names(o)) > 0]
    return(parse_named(o, prefix, abs))
}

parse_snakemake <- function(snakemake){
    # Parse slots from global snakemake object into separate objects
    parse_snakemake_slot(snakemake, "params", "", FALSE)
    parse_snakemake_slot(snakemake, "wildcards", "", FALSE)
    parse_snakemake_slot(snakemake, "input", "inpath", FALSE)
    parse_snakemake_slot(snakemake, "output", "outpath", TRUE)
    parse_snakemake_slot(snakemake, "log", "logpath", FALSE)
}

#-----------------------------------------------------------------------------
# Log script progress
#-----------------------------------------------------------------------------

clear_log <- function(path=get("logpath")){
    # Delete a pre-existing log file
    if (file.exists(path)) file.remove(path)
}

write_log <- function(..., newline=TRUE, path=get("logpath")){
    # Append messages to log file
    cat(c(...), ifelse(newline, "\n", ""), file=path, append=TRUE)
}

log_done <- function(newline=FALSE, path=get("logpath")) {
    write_log("done.", ifelse(newline, "\n", ""), path=path)
}
