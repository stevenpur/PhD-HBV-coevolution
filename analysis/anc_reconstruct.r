library(tidyverse)
library(ape)
library(seqinr)
library(progress)
library(doSNOW)

setwd("~/hbv_covar3")
genotype <- "D"
genes <- c("C", "P", "S", "X")
ncores <- 10

# output file
anc_restruct_wrt_file <- paste0("./analysis/hbv", genotype, "_aa_anc_restruct.rds")

# input file
tree_file <- paste0("./analysis/phylo_build/RAxML_bestTree.HBV", genotype, "_rooted")
gene_msa_files <- map(genes, function(x) {
    paste0("./amino_acid/HBV_gene", x, "_AAseqs_noShift_stopToDot_mafft.fasta")
})
names(gene_msa_files) <- genes

# load data
tree <- read.tree(tree_file)
gene_msa <- map(genes, function(x) {
    do.call(rbind, read.fasta(gene_msa_files[[x]]))
})
names(gene_msa) <- genes

# anc reconstruction for each varying site
message("performing ancestral state reconstruction...")
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
    # add states informatino for tip nodes to anc restruct result
    tip_nodes_lik_anc <- map(site_aas_matched, function(x) {
        ifelse(colnames(ace_result$lik.anc) == x, 1, 0)
    })
    tip_nodes_lik_anc <- do.call(rbind, tip_nodes_lik_anc)
    ace_result$lik.anc <- rbind(tip_nodes_lik_anc, ace_result$lik.anc)
    anc_restruct[[site_id]] <- ace_result
}


saveRDS(anc_restruct, anc_restruct_wrt_file)