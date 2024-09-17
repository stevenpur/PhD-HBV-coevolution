# simulate the tree
library(ape)
library(tidyverse)

# get the population size
pop_sizes <- c(100, 1000, 2000)
run_size <- 100

# simulate the tree using coalescent model
trees <- map(pop_sizes, function(pop_size){
    map(1:run_size, function(run){
        tree <- rcoal(pop_size)
        return(tree)
    })
})
trees <- unlist(trees, recursive = F)

# rescale the trees to the same total branch length
total_branch_lengths <- map_dbl(trees, ~sum(.x$edge.length))
mean_total_branch_length <- mean(total_branch_lengths)
trees_rescaled <- map(trees, function(tree){
    tree$edge.length <- tree$edge.length / sum(tree$edge.length) * mean_total_branch_length
    return(tree)
})

# check the total branch length
total_branch_lengths_rescaled <- map_dbl(trees_rescaled, ~sum(.x$edge.length))
print(total_branch_lengths_rescaled)

# write the scaled and unscaled trees to file
map(1:length(trees), function(i){
    tree <- trees[[i]]
    ind <- i %% run_size + 1
    outfile <- paste0("~/hbv_covar3/analysis/sim_seq/simseq_N", length(tree$tip), "_", ind, ".tree")
    write.tree(tree, outfile)
})


map(1:length(trees_rescaled), function(i){
    tree <- trees_rescaled[[i]]
    ind <- i %% run_size + 1
    outfile <- paste0("~/hbv_covar3/analysis/sim_seq/simseq_N", length(tree$tip), "_", ind, "_rescaled.tree")
    write.tree(tree, outfile)
})
