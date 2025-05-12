# Load required libraries
suppressPackageStartupMessages({
    library(optparse)
    library(tidyverse)
    library(ape)
})

# Parse command line arguments
option_list <- list(
    make_option(c("--pop"), type = "character", help = "Population sizes to simulate (comma-separated)"),
    make_option(c("--runs"), type = "integer", default = 100, help = "Number of runs per population size"),
    make_option(c("--rescale"), action="store_true", default = FALSE, help = "Whether to rescale trees to same total branch length"),
    make_option(c("--no-rescale"), action="store_false", dest="rescale", help = "Disable tree rescaling"),
    make_option(c("--outdir"), type = "character", default = "./", help = "Output directory for simulation results")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Convert population sizes from string to numeric vector
pop_sizes <- as.numeric(strsplit(opt$pop, ",")[[1]])
run_size <- opt$runs
do_rescale <- opt$rescale
outdir <- opt$outdir

# simulate the tree using coalescent model
trees <- map(pop_sizes, function(pop_size){
    map(1:run_size, function(run){
        tree <- rcoal(pop_size)
        return(tree)
    })
})
trees <- unlist(trees, recursive = F)

# rescale the trees to the same total branch length if requested
if (do_rescale) {
    total_branch_lengths <- map_dbl(trees, ~sum(.x$edge.length))
    mean_total_branch_length <- mean(total_branch_lengths)
    trees_rescaled <- map(trees, function(tree){
        tree$edge.length <- tree$edge.length / sum(tree$edge.length) * mean_total_branch_length
        return(tree)
    })
}

# write the trees to file
if (do_rescale) {
    map(1:length(trees_rescaled), function(i){
        tree <- trees_rescaled[[i]]
        ind <- i %% run_size + 1
        outfile <- paste0(outdir, "/simseq_N", length(tree$tip), "rs_", ind, ".tree")
        write.tree(tree, outfile)
    })
} else {
    map(1:length(trees), function(i){
        tree <- trees[[i]]
        ind <- i %% run_size + 1
        outfile <- paste0(outdir, "/simseq_N", length(tree$tip), "_", ind, ".tree")
        write.tree(tree, outfile)
    })
}