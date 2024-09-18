# This script calculates the allele switch counts for each site in a given genotype of the hepatitis B virus (HBV).
# It requires the following libraries: tidyverse, ape, phangorn, foreach, doParallel.

# The script takes command line arguments:
# - genotype: the genotype of the HBV
# - anc_restruct_file: the file containing the ancestral reconstruction data
# - wrt_file: the output file to save the allele switch counts

# The script performs the following steps:
# 1. Sets the working directory to "~/hbv_covar3/"
# 2. Reads the ancestral reconstruction data from the anc_restruct_file
# 3. Reads the phylogenetic tree from the tree_file
# 4. Registers parallel processing with the specified number of cores (ncores)
# 5. Iterates over each site in the ancestral reconstruction data and calculates the allele switch counts
# 6. Saves the allele switch counts to the wrt_file using the saveRDS function

# Usage: Rscript get_allele_homoplasy.r <genotype> <anc_restruct_file> <wrt_file>

# Note: The script assumes that the necessary input files are present in the specified paths.

library(tidyverse)
library(ape)
library(phangorn)
library(foreach)
library(doParallel)

#########################
# Parameters to change  #
#########################
args <- commandArgs(trailingOnly = TRUE)

setwd("~/hbv_covar3/")
genotype <- args[1]
allele_presense_threshold <- 0.7
allele_empty_threshold <- 0.3
ncores <- 10

# output file
wrt_file <- args[3]
# wrt_file <- paste0("./analysis/hbv", genotype, "_allele_switch_cnt.rds")

# input file
anc_restruct_file <- args[2]
# tree_file <- paste0("./analysis/phylo_build/RAxML_bestTree.HBV", genotype, "_rooted")

# load file
anc_restruct <- readRDS(anc_restruct_file)
tree <- read.tree(tree_file)

################
# Run the code #
################
# for each allele, where did the substitution occur?
registerDoParallel(ncores)
allele_switch <- foreach(
    site_id = names(anc_restruct),
    .combine = "c"
) %dopar% {
    message(site_id)
    gene <- gsub("gene", "", strsplit(site_id, "_")[[1]][1])
    site <- as.numeric(gsub("gene._site", "", site_id))
    cur_anc_restruct <- anc_restruct[[site_id]]$lik.anc
    cur_switch_allele_tb <- matrix(NA,
        nrow = nrow(cur_anc_restruct),
        ncol = ncol(cur_anc_restruct)
    )
    colnames(cur_switch_allele_tb) <- colnames(cur_anc_restruct)
    for (node in seq_len(tree$Nnode + Ntip(tree))) {
        # get current node state
        cur_state <- cur_anc_restruct[node, ]
        if (NA %in% cur_state) next
        # get parent node state
        parent_node <- Ancestors(tree, node, "parent")
        if (parent_node == 0) next
        parent_state <- cur_anc_restruct[parent_node, ]
        # identify switch
        is_disappear <- (parent_state > allele_presense_threshold &
            cur_state < allele_empty_threshold)
        is_appear <- (parent_state < allele_empty_threshold &
            cur_state > allele_presense_threshold)
        is_switch <- is_disappear | is_appear
        cur_switch_allele_tb[node, ] <- is_switch
    }
    result <- list()
    result[[site_id]] <- cur_switch_allele_tb
    result
}


saveRDS(allele_switch, wrt_file)