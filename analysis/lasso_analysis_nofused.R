library(tidyverse)
library(ape)
library(seqinr)
library(glmnet)
library(SIS)
library(doParallel)
library(foreach)
library(itertools)

args <- commandArgs(trailingOnly = TRUE)
ncores <- 10
setwd("~/hbv_covar3/")
lg <- list()
lg$genotype <- args[1]
lg$genes <- c("C", "P", "S", "X")

# input files
lg$allele_switch_file <- paste0("./analysis/hbv", lg$genotype, "_allele_switch_cnt.rds")
lg$tree_file <- paste0("./analysis/phylo_build/RAxML_bestTree.HBV", lg$genotype, "_rooted")
lg$ovp_matrix_file <- paste0("./analysis/hbv", lg$genotype, "_overlap.rds")
lg$gene_msa_files <- map(lg$genes, function(x) {
    paste0("./amino_acid/HBV_gene", x, "_AAseqs_noShift_stopToDot_mafft.fasta")
})
names(lg$gene_msa_files) <- lg$genes

# load data
allele_switch <- readRDS(lg$allele_switch_file)
tree <- read.tree(lg$tree_file)
ovp_matrix <- readRDS(lg$ovp_matrix_file)
gene_msa <- map(lg$genes, function(x) {
    do.call(rbind, read.fasta(lg$gene_msa_files[[x]]))
})
names(gene_msa) <- lg$genes

# filter out alleles that do not need to be considered
allele_switch_qced <- allele_switch
rm_ind <- vector()
lg$gap_filter_record <- vector()
lg$switch_cnt_allele_filtered <- vector()
lg$switch_cnt_site_filtered <- vector()
for (i in seq_len(length(allele_switch))) {
    site_id <- names(allele_switch)[i]
    # can't have 99% site as gap
    gene <- gsub("gene", "", strsplit(site_id, "_")[[1]][1])
    site <- as.numeric(gsub("gene._site", "", site_id))
    site_aa <- gene_msa[[gene]][, site]
    site_aa_matched <- site_aa[tree$tip.label]
    site_freq <- table(site_aa) / length(site_aa)
    if ("-" %in% names(site_freq)) {
        if (site_freq["-"] > 0.99) {
            rm_ind <- c(rm_ind, i)
            lg$gap_filter_record <- c(lg$gap_filter_record, site_id)
            next
        }
    }
    # an allele need to have at least 3 switches
    cur_allele_switch <- allele_switch_qced[[site_id]]
    switch_cnt <- apply(cur_allele_switch, 2, function(x) sum(x, na.rm = T))
    if (length(which(switch_cnt <= 2)) != 0) {
        rm_id <- paste0(
            site_id,
            "_",
            colnames(cur_allele_switch)[which(switch_cnt <= 2)]
        )
        lg$switch_cnt_filter_record <- c(lg$switch_cnt_filter_record, rm_id)
        allele_switch_qced[[site_id]] <- cur_allele_switch[, switch_cnt >= 3, drop = F]
    }
    if (ncol(allele_switch_qced[[site_id]]) == 0) {
        lg$switch_cnt_site_filtered <- c(lg$switch_cnt_site_filtered, site_id)
        rm_ind <- c(rm_ind, i)
    }
}
if (length(rm_ind) > 0) {
    allele_switch_qced <- allele_switch_qced[-rm_ind]
}

site_switch <- map(allele_switch_qced, function(x) {
    apply(x, 1, any)
})
nt_msa_len <- 3665

# get tree length
nodes_n <- tree$Nnode + length(tree$tip.label)
root_node <- setdiff(1:nodes_n, tree$edge[, 2])
edge_length <- rep(NA, nodes_n)
edge_length[tree$edge[, 2]] <- tree$edge.length

warning_log <- vector()
registerDoParallel(ncores)
iterate_chunks <- isplitVector(seq_len(length(site_switch)), chunks = ncores)

result <- foreach(chunk = iterate_chunks) %dopar% {
    cur_result <- list()
    for (i in chunk) {
        site_id <- names(site_switch)[i]
        message(i)
        ovp_sites <- c(colnames(ovp_matrix)[which(ovp_matrix[site_id, ] != 0)], site_id)
        # fuse the nodes together
        site_switch_fused <- site_switch
        dep_sites <- names(site_switch_fused)[-which(names(site_switch_fused) %in% ovp_sites)]
        dep_sites <- do.call(cbind, site_switch_fused[dep_sites])
        dep_sites <- apply(dep_sites, 2, as.numeric)
        data <- cbind(site_switch[[site_id]], edge_length, dep_sites)
        # remove cols/rows with too many NA
        rm_row <- unlist(apply(data, 1, function(x) {
            freq <- sum(is.na(x)) / length(x)
            if (freq > 0.1) {
                return(T)
            } else {
                return(F)
            }
        }))

        rm_col <- unlist(apply(data, 2, function(x) {
            freq <- sum(is.na(x)) / length(x)
            if (freq > 0.1) {
                return(T)
            } else {
                return(F)
            }
        }))
        if (length(rm_row)) {
            data <- data[!(rm_row), !(rm_col)]
        }
        y <- data[, 1]
        x <- data[, -1]
        mysd <- function(y) sqrt(sum((y - mean(y))^2) / length(y))
        sx <- scale(x, scale = apply(x, 2, mysd))
        start_lambda <- max(abs(colSums(sx * y))) / nrow(x)
        # add 100 lambdas before start lambda for warm up
        lambdas <- map_dbl(1:100, function(i) start_lambda * (1 / 0.99)^(i - 1))
        lambdas <- c(lambdas, map_dbl(1:687, function(i) start_lambda * (0.99)^(i - 1)))
        lambdas <- sort(lambdas)
        cur_result[[site_id]] <- glmnet(data[, -1], data[, 1], family = "binomial", lambda = lambdas)
        # select beta via ebic
        ebic <- map_dbl(1:length(cur_result[[site_id]]$lambda), function(j) {
            lambda <- cur_result[[site_id]]$lambda[j]
            cur_ebic <- deviance(cur_result[[site_id]])[j] +
                sum(cur_result[[site_id]]$beta[, j] != 0) * log(nrow(data)) +
                2 * log(choose(ncol(data[, -1]), sum(cur_result[[site_id]]$beta[, j] != 0)))
        })
        min.ind <- which.min(ebic)
        beta_col <- cur_result[[site_id]]$beta[, min.ind]
        beta <- names(beta_col)[which(beta_col != 0)]
        message(beta)
        cur_result[[site_id]]$selected_beta <- beta
    }
    return(cur_result)
}

result <- unlist(result, recursive = F)


saveRDS(result, paste0("./analysis/lasso_covar_hbv", lg$genotype, ".rds"))