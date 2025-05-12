library(tidyverse)
library(reshape2)
library(gridExtra)
library(ggplot2)
setwd("~/hbv_covar3/analysis/sim_seq/")

#%% define functions
get_lnfu_comb <- function(files, is_null = FALSE){
    # get the running parameter for each result
    files <- gsub(".*/", "", files)
    files_strat <- map_chr(strsplit(files, "_"), 2)
    files_strat <- strsplit(files_strat, '[lnfu]')
    strat_tb <- do.call(rbind, files_strat)[, -1]
    strat_tb <- as.data.frame(apply(strat_tb, 2, as.numeric))
    if (is_null) {
        colnames(strat_tb) <- c("l", "n", "u")
        l_unique <- unique(strat_tb$l)
        n_unique <- unique(strat_tb$n)
        u_unique <- unique(strat_tb$u)
        lnfu_comb <- expand.grid(l_unique, n_unique, u_unique)
        colnames(lnfu_comb) <- c("l", "n", "u")
        lnfu_comb$param_set <- paste0("l", lnfu_comb$l, "n", lnfu_comb$n, "u", lnfu_comb$u)
    } else {
        colnames(strat_tb) <- c("l", "n", "f", "u")
        l_unique <- unique(strat_tb$l)
        n_unique <- unique(strat_tb$n)
        f_unique <- unique(strat_tb$f)
        u_unique <- unique(strat_tb$u)
        lnfu_comb <- expand.grid(l_unique, n_unique, f_unique, u_unique)
        colnames(lnfu_comb) <- c("l", "n", "f", "u")
        lnfu_comb$param_set <- paste0("l", lnfu_comb$l, "n", lnfu_comb$n, "f", lnfu_comb$f, "u", lnfu_comb$u)
    }
    return(lnfu_comb)
}

get_tp_fp <- function(tb, coev_pair, indir_coev_pair, l){
    if (class(tb) == "logical"){
        tp_cnt <- 0
        fp1_cnt <- 0
        fp2_cnt <- 0
    } else {
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
    }
    fn_cnt <- length(coev_pair) - tp_cnt
    tn_cnt <- choose(l, 2) - tp_cnt - fp1_cnt - fp2_cnt - fn_cnt
    result <- c(tp_cnt, fp1_cnt, fp2_cnt, fn_cnt, tn_cnt)
    names(result) <- c("tp", "fp1", "fp2", "fn", "tn")
    return(result)
}

get_metrics <- function(cm) {
    cm$fp <- sum(cm$fp1, cm$fp2)
    accuracy <- (cm$tp + cm$tn) / sum(unlist(cm))
    recall <- cm$tp / (cm$tp + cm$fn)
    precision <- cm$tp / (cm$tp + cm$fp)
    specificity <- cm$tn / (cm$tn + cm$fp)
    f1 <- 2 * precision * recall / (precision + recall)
    balanced_accuracy <- (recall + specificity) / 2
    dir_indir_ratio <- cm$fp1 / cm$tp
    result <- c(accuracy, recall, precision, specificity, balanced_accuracy, balanced_accuracy, dir_indir_ratio)
    names(result) <- c("accuracy", "recall", "precision", "specificity", "f1", "balanced_accuracy", "dir_indir_ratio")
    return(result)
}


#%% result files for treeso
pattern <- 'l100'
folders <- list.dirs(full.names = FALSE, recursive = FALSE)
param_set <- folders[str_detect(folders, pattern)]

pattern <- "simresult.*txt"
result_files <- unlist(map(param_set, ~list.files(.x, pattern = pattern, full.names = TRUE)))

if (length(result_files) == 0) {
  stop("No result files found")
}
coev_pair <- c("site1:site2", "site2:site3") # the true coevolving pairs
indir_coev_pair <- c("site1:site3") # the indirect coevolution
true_dir_assoc_cnt <- 2
true_indir_assoc_cnt <- 1

#%% get the performance of treeso for each parameter set 

lnf_comb <- get_lnfu_comb(result_files)
result <- list()

for (i in 1:nrow(lnf_comb)) {
    param_set <- lnf_comb$param_set[i]
    result[[param_set]] <- list()
    message(paste0('calculating for parameter set: ', param_set))
    param_set_pattern <- paste0(param_set, "_")
    cur_files <- result_files[grepl(param_set_pattern, result_files)]
    message(paste0('number of files: ', length(cur_files)))
    tp_fp <- map(cur_files, function(file) {
        # if there is no line, then there is no TP nor FP
        if (file.info(file)$size < 2) {
            cur_tp_fp <- get_tp_fp(NA, coev_pair, indir_coev_pair, lnf_comb$l[i])
        } else {
            # get the pair id from the result files
            tb <- read.table(file, sep = " ", header = F)
            colnames(tb) <- c("siteA", "siteB")
            tb$pair_id <- map_chr(1:nrow(tb), function(x) {
                paste(sort(c(tb[x, 1], tb[x, 2])), collapse = ":")
            })
            # use the pair_id to find the TP and FP
            cur_tp_fp <- get_tp_fp(tb, coev_pair, indir_coev_pair, lnf_comb$l[i])
        }
        return(cur_tp_fp)
    })
    # column-wise sum of the TP and FP
    tp_fp_sum <- as.list(colSums(do.call(rbind, tp_fp)))
    result[[param_set]] <- get_metrics(tp_fp_sum)
}

roc_result <- as.data.frame(do.call(rbind, result))

# roc_result <- read.table("./treeso_roc_result.txt", sep = "\t", header = T)
write.table(roc_result, './treeso_roc_result.txt', sep = '\t')

#%% investigate more into why treeso is failing
lnf_comb <- "l100n1000f100u0.01" 
cur_files <- result_files[grepl(lnf_comb, result_files)]
for (i in 1:100){
    file_name <- paste0("test//simresult_", lnf_comb, "_", i, ".txt")
    if (! (file_name %in% cur_files)) {
        message(paste0("file ", file_name, " is not in the list"))
    } 
}




#%% save the results
setwd("~/hbv_covar3/analysis/sim_seq/")
write.table(roc_result, file = "roc_result.txt", sep = "\t", row.names = T)

#%% read the results of the dca results
pattern <- 'l100'
folders <- list.dirs(full.names = FALSE, recursive = FALSE)
param_set <- folders[str_detect(folders, pattern)]

pattern <- "simresult.*dca.csv"
dca_files <- unlist(map(param_set, ~list.files(.x, pattern = pattern, full.names = TRUE)))

dca_files1 <- dca_files[which(str_detect(dca_files, "f.*mfdca.*"))]
dca_files2 <- setdiff(dca_files, dca_files1)

lnf_comb1 <- get_lnfu_comb(dca_files1, is_null = FALSE)
lnf_comb2 <- get_lnfu_comb(dca_files2, is_null = TRUE)

lnf_comb

#%% Getting the results of DCA
#%% Method 1: get the top L pairs
l = 100 
result <- list()
for (i in 1:nrow(lnf_comb1)){
    param_set <- lnf_comb1$param_set[i]
    message(paste0('calculating for parameter set: ', param_set))
    cur_files <- dca_files1[grepl(param_set, dca_files1)]
    tp_fp <- map(cur_files, function(file){
        if (file.exists(file) == F | file.info(file)$size < 2) {
            cur_tp_fp <- get_tp_fp(NA, coev_pair, indir_coev_pair, lnf_comb1$l[i])
        } else {
            tb <- read.csv(file, header = F)
            colnames(tb) <- c("siteA", "siteB", "DC_strength")
            tb <- tb[1:2, ]
            tb$siteA <- paste0("site", tb$siteA+1)
            tb$siteB <- paste0("site", tb$siteB+1)
            tb$pair_id <- map_chr(1:nrow(tb), function(x) {
                paste(sort(c(tb[x, 1], tb[x, 2])), collapse = ":")
            })
            if (length(tb) == 0) {
                cur_tp_fp <- get_tp_fp(NA, coev_pair, indir_coev_pair, lnf_comb1$l[i])
            } else {
                cur_tp_fp <- get_tp_fp(tb, coev_pair, indir_coev_pair, lnf_comb1$l[i])
            }
        }
        return(cur_tp_fp)
    })
    tp_fp_sum <- as.list(colSums(do.call(rbind, tp_fp)))
    result[[param_set]] <- get_metrics(tp_fp_sum)
}

dca_result <- as.data.frame(do.call(rbind, result))

#%% Method 2: get the performance of dca for each parameter set
# first we want to get the null distribution
null_thresh <- list()
coeffs <- list()
for (i in 1:nrow(lnf_comb2)) {
    param_set <- lnf_comb2$param_set[i]
    message(paste0('calculating for parameter set: ', param_set))
    null_thresh[[param_set]] <- vector() 
    cur_files <- dca_files2[grepl(param_set, dca_files2)]
    coeffs[[param_set]] <- map(cur_files, function(file){
        if (file.exists(file) == F) {

        } else if (file.info(file)$size < 2) {
            return(c())
        } else {
            tb <- read.csv(file, header = F)
            return(tb[, 3])
        }
    })
    coeffs[[param_set]] <- unlist(coeffs[[param_set]])
    coeffs[[param_set]] <- sort(coeffs[[param_set]])
    thresh_ind <- round(length(coeffs[[param_set]]) * 0.95)
    null_thresh[[param_set]] <- coeffs[[param_set]][thresh_ind]
}

dca_result <- list()
for (i in 1:nrow(lnf_comb1)) {
    param_set <- lnf_comb1$param_set[i]
    null_param_set <- gsub("f[0-9.]+", "", param_set)
    message(paste0('calculating for parameter set: ', param_set))
    dca_result[[param_set]] <- list()
    cur_files <- dca_files[grepl(param_set, dca_files)]
    tp_fp <- map(cur_files, function(file){
        if (file.info(file)$size < 2) {
            cur_tp_fp <- get_tp_fp(NA, coev_pair, indir_coev_pair, lnf_comb1$l[i])
        } else {
            tb <- read.csv(file, header = F)
            colnames(tb) <- c("siteA", "siteB", "DC_strength")
            tb <- tb[tb$DC_strength >= null_thresh[[null_param_set]], ]
            tb$siteA <- paste0("site", tb$siteA+1)
            tb$siteB <- paste0("site", tb$siteB+1)
            tb$pair_id <- map_chr(1:nrow(tb), function(x) {
                paste(sort(c(tb[x, 1], tb[x, 2])), collapse = ":")
            })
            if (length(tb) == 0) {
                cur_tp_fp <- get_tp_fp(NA, coev_pair, indir_coev_pair, lnf_comb1$l[i])
            } else {
                cur_tp_fp <- get_tp_fp(tb, coev_pair, indir_coev_pair, lnf_comb1$l[i])
            }
        }
            return(cur_tp_fp)
        })
    tp_fp_sum <- as.list(colSums(do.call(rbind, tp_fp)))
    dca_result[[param_set]] <- get_metrics(tp_fp_sum)
}

dca_roc_result <- as.data.frame(do.call(rbind, dca_result))

#%% save the results
setwd("/users/bag/hlq763/hbv_covar3/analysis/sim_seq/")
write.table(dca_roc_result, file = "dca_roc_result.txt", sep = "\t", row.names = T)

#%% match the two tables
roc_result$param_set <- rownames(roc_result)
dca_result$param_set <- rownames(dca_result)

combined_roc <- merge(roc_result, dca_result, by = "param_set")

write.table(combined_roc, file = "combined_roc_result.temp2.txt", sep = "\t", row.names = F)


