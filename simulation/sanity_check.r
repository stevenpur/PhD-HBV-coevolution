#%% import libraries
library(tidyverse)

#%% check to see if all simulation result is available 
sim_dir <- '~/hbv_covar3/analysis/sim_seq/' 
setwd(sim_dir)

u <- c(0.01, 0.1, 1, 10, 100)
f <- c(2, 10, 100)
n <- c(100, 1000, 2000)
run_size <- 100

param_comb <- expand.grid(u = u, f = f, n = n)
param_set <- map_chr(1:nrow(param_comb), function(i) {
    paste0('l100', 'n', param_comb[i, 'n'], 'f', param_comb[i, 'f'], 'u', param_comb[i, 'u'])
})

for (param in param_set) {
    for (ind in 1:run_size){
        simseq_file <- paste0(param, '/simseq_', param, '_', ind, '.fasta')
        if (!file.exists(simseq_file)) {
            print(paste0('Missing: ', simseq_file))
        }
    }
}

for (param in param_set) {
    for (ind in 1:run_size){
        simseq_file <- paste0(param, '/simseq_', param, '_', ind, 'rescaled.fasta')
        if (!file.exists(simseq_file)) {
            print(paste0('Missing: ', simseq_file))
        }
    }
}

simresult_files <- list.files(sim_dir, pattern = 'simresult_')

param_set <- map_chr(simresult_files, ~str_match(.x, 'simresult_(.*)\\_.*.txt')[,2])

param_set_unique <- unique(param_set)
param_set_cnt <- vector('list', length(param_set_unique))

for (i in 1:length(param_set_unique)) {
    cur_param_set <- param_set_unique[i]
    pattern = paste0('_', cur_param_set, '_')
    cur_sim_files <- simresult_files[str_detect(simresult_files, pattern)]
    param_set_cnt[[i]] <- c(cur_param_set, as.character(length(cur_sim_files)))
}

#%% for the missing result, check the log files

missing_si 
