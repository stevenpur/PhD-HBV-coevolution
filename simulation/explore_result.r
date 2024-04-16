library(tidyverse)
library(reshape2)
library(gridExtra)
library(ggplot2)
setwd("~/hbv_covar3/analysis/sim_seq/")

#%% define functions
get_lnfu_comb <- function(files){
    # get the running parameter for each result
    files <- gsub(".*/", "", files)
    files_strat <- map_chr(strsplit(files, "_"), 2)
    files_strat <- strsplit(files_strat, '[lnfu]')
    strat_tb <- do.call(rbind, files_strat)[, -1]
    strat_tb <- as.data.frame(apply(strat_tb, 2, as.numeric))
    colnames(strat_tb) <- c("l", "n", "f", "u")
    l_unique <- unique(strat_tb$l)
    n_unique <- unique(strat_tb$n)
    f_unique <- unique(strat_tb$f)
    u_unique <- unique(strat_tb$u)
    lnfu_comb <- expand.grid(l_unique, n_unique, f_unique, u_unique)
    colnames(lnfu_comb) <- c("l", "n", "f", "u")
    lnfu_comb$param_set <- paste0("l", lnfu_comb$l, "n", lnfu_comb$n, "f", lnfu_comb$f, "u", lnfu_comb$u)
    return(lnfu_comb)
}

get_tp_fp <- function(tb, coev_pair, indir_coev_pair){
    if (!'pair_id' %in% colnames(tb)) {
        stop("pair_id not in the tb")
    }
    pair_id_unique <- unique(tb$pair_id)
    el_ind <- grep("edge_length", pair_id_unique)
    # we are not interested in edge_length signal
    if (length(el_ind) > 0) {
        pair_id_unique <- pair_id_unique[-el_ind]
    }
    tp_cnt <- sum(pair_id_unique %in% coev_pair)
    fp1_cnt <- sum(pair_id_unique %in% indir_coev_pair)
    fp2_cnt <- length(pair_id_unique) - tp_cnt - fp1_cnt
    return (c(tp_cnt, fp1_cnt, fp2_cnt))
}

get_tpr_fpr <- function(true_dir_assoc_cnt, true_indir_assoc_cnt, l){
    true_no_assoc_cnt <- choose(l, 2)
    avg_tpr <- sum(tp_fp$coev_found)/(true_dir_assoc_cnt*nrow(tp_fp))
    avg_fpr <- (sum(tp_fp$indir_coev_found) + sum(tp_fp$false_coev_found))/(true_no_assoc_cnt * nrow(tp_fp))
    avg_indir <- sum(tp_fp$indir_coev_found) / nrow(tp_fp)
    return(c(avg_tpr, avg_fpr, avg_indir))
}


#%% result files
pattern = "simresult.*txt"
result_files <- list.files("./", pattern = pattern)

if (length(result_files) == 0) {
  stop("No result files found")
}
coev_pair <- c("site1:site2", "site2:site3") # the true coevolving pairs
indir_coev_pair <- c("site1:site3") # the indirect coevolution
true_dir_assoc_cnt <- 2
true_indir_assoc_cnt <- 1

#%% get the count of each files
lnf_comb <- get_lnfu_comb(result_files)
result <- list()

for (i in 1:nrow(lnf_comb)) { 
    param_set <- lnf_comb$param_set[i]
    result[[param_set]] <- list()
    message(paste0('calculating for parameter set: ', param_set))
    cur_files <- result_files[grepl(param_set, result_files)]
    tp_fp <- map(cur_files, function(file) {
        # if there is no line, then there is no TP nor FP
        if (file.info(file)$size < 2) {
            return(c(0, 0, 0))
        } else {
            # get the pair id from the result files
            tb <- read.table(file, sep = " ", header = F)
            colnames(tb) <- c("siteA", "siteB")
            tb$pair_id <- map_chr(1:nrow(tb), function(x) {
                paste(sort(c(tb[x, 1], tb[x, 2])), collapse = ":")
            })
            # use the pair_id to find the TP and FP
            tp_fp <- get_tp_fp(tb, coev_pair, indir_coev_pair)
        }
        return(tp_fp)
    })
    tp_fp <- as.data.frame(do.call(rbind, tp_fp))
    colnames(tp_fp) <- c("coev_found", "indir_coev_found", "false_coev_found")
    result[[param_set]] <- get_tpr_fpr(true_dir_assoc_cnt, true_indir_assoc_cnt, lnf_comb$l[i])
}

roc_result <- as.data.frame(do.call(rbind, result)
colnames(roc_result) <- c("avg_tpr", "avg_fpr", "avg_indir")

#%% plot the results in contour plot
# Replace 'path_to_your_data.csv' with the actual path to your data file
data <- roc_result
# Data preparation





# read the results of the dca results
dca_files = list.files("/users/bag/hlq763/hbv_covar3/analysis/sim_seq/archive/", pattern = "simseq_.*mfdca.csv", full.names = T)
dca_result <- list()
lnf_comb2 <- get_lnfu_comb(dca_files)

dca_result <- list()
for (i in 1:nrow(lnf_comb2)) {
    param_set <- lnf_comb2$param_set[i]
    message(paste0('calculating for parameter set: ', param_set))
    dca_result[[param_set]] <- list()
    cur_files <- dca_files[grepl(param_set, dca_files)]
    tp_fp <- map(cur_files, function(file){
        if (file.info(file)$size < 2) {
            return(c(0, 0, 0))
        } else {
            tb <- read.csv(file, header = F)
            colnames(tb) <- c("siteA", "siteB", "DC_strength")
            tb$siteA <- paste0("site", tb$siteA+1)
            tb$siteB <- paste0("site", tb$siteB+1)
            tb$pair_id <- map_chr(1:nrow(tb), function(x) {
                paste(sort(c(tb[x, 1], tb[x, 2])), collapse = ":")
            })
            tp_fp <- get_tp_fp(tb, coev_pair, indir_coev_pair)
        }
            return(tp_fp)
        })
    tp_fp <- as.data.frame(do.call(rbind, tp_fp))
    colnames(tp_fp) <- c("coev_found", "indir_coev_found", "false_coev_found")
    result[[param_set]] <- get_tpr_fpr(true_dir_assoc_cnt, true_indir_assoc_cnt, lnf_comb2$l[i])
}

dca_roc_result <- do.call(rbind, result)
colnames(dca_roc_result) <- c("avg_tpr", "avg_fpr", "avg_indir")


