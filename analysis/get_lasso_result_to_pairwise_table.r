library(tidyverse)
setwd("~/hbv_covar3")

genotype <- "D"
assoc_file <- paste0("./analysis/hbv", genotype, "_glm_edge_test_extraInfo.tsv")
lasso_file <- paste0("./analysis/lasso_covar_hbv", genotype, ".rds")
bs_file <- paste0("./analysis/lasso_covar_hbv_bootstrap", genotype, ".rds")

lasso_result <- readRDS(lasso_file)
glm_assoc <- read.table(assoc_file, header = T, stringsAsFactors = F, sep = "\t")
bs <- readRDS(bs_file)

chosen_betas <- map(lasso_result, function(x) x$selected_beta)
lasso_pair <- unlist(map(1:length(chosen_betas), function(i) {
    if (length(chosen_betas[[i]]) == 0) {
        return(NULL)
    }
    siteAs <- strsplit(names(chosen_betas)[i], ":")[[1]]
    siteBs <- unlist(strsplit(chosen_betas[[i]], ":"))
    cur_pairs_t <- expand.grid(siteAs, siteBs, stringsAsFactors = F)
    cur_pairs <- map_chr(1:nrow(cur_pairs_t), function(j) {
        siteA <- cur_pairs_t[j, 1]
        siteB <- cur_pairs_t[j, 2]
        return(paste(sort(c(siteA, siteB)), collapse = "_"))
    })
    return(cur_pairs)
}))


lasso_pair_freq <- table(lasso_pair)
strong_lasso_pair <- names(lasso_pair_freq)[which(lasso_pair_freq == 2)]
weak_lasso_pair <- names(lasso_pair_freq)[which(lasso_pair_freq == 1)]


assoc_lasso <- glm_assoc
assoc_lasso$adj_pval <- 1
assoc_lasso$adj_pval[which(glm_assoc$site_pair_id %in% strong_lasso_pair)] <- 0.01
assoc_lasso$adj_pval[which(glm_assoc$site_pair_id %in% weak_lasso_pair)] <- 0.001

# get bootstrap value for the associations
assoc_lasso$bs1 <- NA
assoc_lasso$bs2 <- NA
sig_row <- which(assoc_lasso$adj_pval < 0.05)
for (row in sig_row) {
    siteA <- assoc_lasso$siteA[row]
    siteB <- assoc_lasso$siteB[row]
    bs1_cnt <- map_dbl(bs[[siteA]]$selected_beta, function(cur_beta) {
        if (siteB %in% cur_beta) {
            return(1)
        } else {
            return(0)
        }
    })
    assoc_lasso$bs1[row] <- sum(bs1_cnt) / length(bs1_cnt)
    bs2_cnt <- map_dbl(bs[[siteB]]$selected_beta, function(cur_beta) {
        if (siteA %in% cur_beta) {
            return(1)
        } else {
            return(0)
        }
    })
    assoc_lasso$bs2[row] <- sum(bs2_cnt) / length(bs2_cnt)
}


write.table(assoc_lasso, paste0("./analysis/hbv", genotype, "_lasso_nofuse_pair_extraInfo_mspval_hmp.tsv"))