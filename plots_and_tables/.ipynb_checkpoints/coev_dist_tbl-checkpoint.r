library(tidyverse)
library(seqinr)
library(scales)
setwd("~/hbv_covar3")
genotypes <- c("A", "B", "C", "D")
genes <- c("C", "P", "S", "X")

#%% load data
assoc_files <- paste0("./analysis/hbv", genotypes, "_lasso_nofuse_pair_extraInfo_mspval_hmp_domain.tsv")
names(assoc_files) <- genotypes

assocs <- map(assoc_files, function(x) {
    read.table(x, header = T, stringsAsFactors = F)
})
names(assocs) <- genotypes

gene_msa_files <- map(genes, function(x) {
    paste0("./amino_acid/HBV_gene", x, "_AAseqs_noShift_stopToDot_mafft.fasta")
})
names(gene_msa_files) <- genes

gene_msa <- map(genes, function(x) {
    do.call(rbind, read.fasta(gene_msa_files[[x]]))
})
names(gene_msa) <- genes

tree_files <- paste0("./analysis/phylo_build/RAxML_bestTree.HBV", genotypes, "_rooted")
names(tree_files) <- genotypes

message("loading tree files...")
trees <- map(tree_files, read.tree)
names(trees) <- genotypes


domain_pos_file <- "./analysis/protein_domain_pos_aa.rds"
domain_pos <- readRDS(domain_pos_file)

#%% get No. sites, No. variable sites, No. covariation
result <- list()
result_simple <- list()

for (genotype in genotypes) {
    tree <- trees[[genotype]]
    cur_dom_pos <- unlist(domain_pos[[genotype]], recursive = F)
    dom_site_cnt <- list()
    for (i in 1:length(cur_dom_pos)) {
        gene_dom <- names(cur_dom_pos)[i]
        gene <- strsplit(gene_dom, "\\.")[[1]][1]
        dom <- strsplit(gene_dom, "\\.")[[1]][2]
        dom_start <- cur_dom_pos[[i]][1]
        dom_end <- cur_dom_pos[[i]][2]
        samps <- intersect(tree$tip.label, rownames(gene_msa[[gene]]))
        cur_msa <- gene_msa[[gene]][samps, dom_start:dom_end]
        rm_col_cnt <- 0
        for (j in 1:ncol(cur_msa)) {
            if (all(cur_msa[, j] == "-")){
                rm_col_cnt <- rm_col_cnt + 1
            }
        }
        dom_site_cnt[[dom]] <- ncol(cur_msa) - rm_col_cnt   
    }

    assoc <- assocs[[genotype]]
    site_ids <- c(assoc$siteA, assoc$siteB)
    domains <- c(assoc$domain1, assoc$domain2)
    domains_nodup <- domains[which(!duplicated(site_ids))]
    dom_var_tb <- table(domains_nodup)


    assoc_sig <- assoc[which(assoc$adj_pval < 0.05 & 
                            (assoc$bs1 > 0.9 | assoc$bs2 > 0.9)), ]
    site_ids <- c(assoc_sig$siteA, assoc_sig$siteB)
    domains <- c(assoc_sig$domain1, assoc_sig$domain2)
    domains_nodup <- domains[which(!duplicated(site_ids))]
    dom_assoc_tb <- table(domains_nodup)

    domains <- map_chr(strsplit(names(cur_dom_pos), "\\."), 2)

    result[[genotype]] <- vector(length = length(domains))
    result_simple[[genotype]] <- vector(length = length(domains))
    for (i in 1:length(domains)) {
        domain <- domains[i]
        dom_var_cnt <- ifelse(is.na(dom_var_tb[domain]), 0, dom_var_tb[domain])
        dom_assoc_cnt <- ifelse(is.na(dom_assoc_tb[domain]), 0, dom_assoc_tb[domain])
        percent1 <- label_percent(accuracy = 0.1)(dom_var_cnt/dom_site_cnt[[domain]])
        percent2 <- label_percent(accuracy = 0.1)(dom_assoc_cnt/dom_site_cnt[[domain]])
        result[[genotype]][i] <- paste(c(dom_var_cnt, dom_site_cnt[[domain]]), collapse = "/")
        result[[genotype]][i] <- paste0(result[[genotype]][i], 
                                        "(", percent1, ")")
        result[[genotype]][i] <- paste0(result[[genotype]][i],
                                        paste(c(dom_assoc_cnt, dom_site_cnt[[domain]]), collapse = "/"))
        result[[genotype]][i] <- paste0(result[[genotype]][i], 
                                        "(", percent2, ")")

        result_simple[[genotype]][i] <- paste(c(dom_assoc_cnt, dom_site_cnt[[domain]]), collapse = "/")
    }
}

result <- do.call(cbind, result)
rownames(result) <- domains

result_simple <- do.call(cbind, result_simple)
rownames(result_simple) <- domains

write.table(result, "./plots/hbv_domain_covar.txt", sep = "\t", quote = F)
write.table(result_simple, "./plots/hbv_domain_covar_simple.txt", sep = "\t", quote = F)
