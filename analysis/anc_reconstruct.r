# This script performs ancestral state reconstruction for a specified genotype of the hepatitis B virus (HBV).
# It uses the `tidyverse`, `ape`, `seqinr`, `progress`, and `doSNOW` libraries.

#########################
# Parameters to change  #
#########################

# Set the working directory to "~/hbv_covar3"
setwd("~/hbv_covar3")

# Set the number of cores to use for parallel processing
ncores <- 10

# Specify the msa files to read in from
genes <- c("C", "P", "S", "X")
gene_msa <- map(genes, function(x) {
    file <- paste0("./amino_acid/HBV_gene", x, "_AAseqs_noShift_stopToDot_mafft.fasta")
    do.call(rbind, read.fasta(file))
})
names(gene_msa) <- genes

# Specify the output file path for the ancestral state reconstruction results
genotype <- "D"
anc_restruct_wrt_file <- paste0("./analysis/hbv", genotype, "_aa_anc_restruct.rds")

# Specify the tree file
tree_file <- paste0("./analysis/phylo_build/RAxML_bestTree.HBV", genotype, "_rooted")
tree <- read.tree(tree_file)

################
# Run the code #
################

# Perform ancestral state reconstruction for each varying site
site_ids <- unlist(map(genes, function(gene) {
    cur_msa <- gene_msa[[gene]]
    paste0("gene", gene, "_site", seq_len(ncol(cur_msa)))
}))
anc_restruct <- list()
pb <- progress_bar$new(total = length(site_ids))
for (site_id in site_ids) {
    pb$tick() # update progress
    Sys.sleep(1 / length(site_ids))
    gene <- gsub("gene", "", strsplit(site_id, "_")[[1]][1])
    site <- as.numeric(gsub("gene._site", "", site_id))
    site_aas <- gene_msa[[gene]][, site]
    site_aas_matched <- site_aas[tree$tip.label]
    if (length(table(site_aas_matched)) == 1) next
    ace_result <- ace(site_aas_matched,
        tree,
        type = "discrete",
        method = "ML"
    )
    # Add states information for tip nodes to ancestral reconstruction result
    tip_nodes_lik_anc <- map(site_aas_matched, function(x) {
        ifelse(colnames(ace_result$lik.anc) == x, 1, 0)
    })
    tip_nodes_lik_anc <- do.call(rbind, tip_nodes_lik_anc)
    ace_result$lik.anc <- rbind(tip_nodes_lik_anc, ace_result$lik.anc)
    anc_restruct[[site_id]] <- ace_result
}

# Save the ancestral reconstruction results as an RDS file
saveRDS(anc_restruct, anc_restruct_wrt_file)