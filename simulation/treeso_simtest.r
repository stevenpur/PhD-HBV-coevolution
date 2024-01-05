library(ape)
library(phangorn)
library(tidyverse)
library(progress)
library(foreach)
library(doParallel)
library(glmnet)
library(SIS)
library(itertools)

setwd("~/hbv_covar3/analysis/sim_seq/")
args <- commandArgs(trailingOnly = TRUE)
run_id <- args[1]
sim_file <- paste0("./simseq_", run_id, ".txt")
tree_file <- paste0("./simseq_", run_id, ".tree")
outfile <- paste0("./simresult_", run_id, ".txt")

msa <- seqinr::read.fasta(sim_file, seqtype = "DNA")
msa <- do.call(rbind, msa)
tree <- read.tree(tree_file)
colnames(msa) <- paste0("site", 1:ncol(msa))
site_ids <- colnames(msa)

##################################
# Ancestral state reconstruction #
##################################

anc_restruct <- list()
pb <- progress_bar$new(total = length(site_ids))
for (site_id in site_ids) {
    pb$tick() # update progress
    Sys.sleep(1 / length(site_ids))
    site <- as.numeric(gsub("site", "", site_id))
    site_aas <- msa[, site]
    names(site_aas) <- rownames(msa)
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

#########################
# Record mutation event #
#########################
allele_presense_threshold <- 0.7
allele_empty_threshold <- 0.3
ncores <- 10

registerDoParallel(ncores)
allele_switch <- foreach(
    site_id = names(anc_restruct),
    .combine = "c"
) %dopar% {
    message(site_id)
    site <- as.numeric(gsub("site", "", site_id))
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


###########################
# Lasso regression search #
###########################
lg <- list()
ncores <- 10

# filter out alleles that do not need to be considered
allele_switch_qced <- allele_switch
rm_ind <- vector()
lg$gap_filter_record <- vector()
lg$switch_cnt_allele_filtered <- vector()
lg$switch_cnt_site_filtered <- vector()
for (i in seq_len(length(allele_switch))) {
    site_id <- names(allele_switch)[i]
    # can't have 99% site as gap
    site <- as.numeric(gsub("site", "", site_id))
    site_aa <- msa[, site]
    names(site_aa) <- rownames(msa)
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

# get tree length
nodes_n <- tree$Nnode + length(tree$tip.label)
root_node <- setdiff(1:nodes_n, tree$edge[, 2])
edge_length <- rep(NA, nodes_n)
edge_length[tree$edge[, 2]] <- tree$edge.length
registerDoParallel(ncores)
iterate_chunks <- isplitVector(seq_len(length(site_switch)), chunks = ncores)

result <- foreach(chunk = iterate_chunks) %dopar% {
    cur_result <- list()
    for (i in chunk) {
        message(i)
        site_id <- names(site_switch)[i]
        # fuse the nodes together
        site_switch_fused <- site_switch
        dep_sites <- names(site_switch_fused)[-which(names(site_switch_fused) %in% site_id)]
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
        x <- data[, -1, drop = F]
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
        # message(paste0("beta_col: ", beta_col))
        beta <- names(beta_col)[which(beta_col != 0)]
        cur_result[[site_id]]$selected_beta <- beta
    }
    return(cur_result)
}

result <- unlist(result, recursive = F)


# check result
assoc_pair <- list()
lst_ind <- 1
for (i in 1:length(result)) {
    cur_assoc_sites <- result[[i]]$selected_beta
    if (length(cur_assoc_sites) != 0) {
        for (cur_assoc_site in cur_assoc_sites) {
            assoc_pair[[lst_ind]] <- c(names(result)[i], cur_assoc_site)
            lst_ind <- lst_ind + 1
        }
    }
}
assoc_pair <- do.call(rbind, assoc_pair)
write.table(assoc_pair, outfile, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
