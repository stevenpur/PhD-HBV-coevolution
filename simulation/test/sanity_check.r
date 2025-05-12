#%% import libraries
library(tidyverse)

#%% check to see if all simulation result is available 
sim_dir <- '~/hbv_covar3/analysis/sim_seq/' 
setwd(sim_dir)

# read the files the directory
pattern = 'l100'
folders <- list.dirs(full.names = FALSE, recursive = FALSE)
folders <- folders[str_detect(folders, pattern)]

# get the parameter set
folders_param <- as.data.frame(do.call(rbind, map(folders, ~str_match(.x, 'l100n(.*)f(.*)u(.*)')[,2:4])))
colnames(folders_param) <- c('n', 'f', 'u')

u <- as.numeric(unique(folders_param$u))
f <- as.numeric(unique(folders_param$f))
n <- as.numeric(unique(folders_param$n))
run_size <- 100

param_comb <- expand.grid(u = u, f = f, n = n)
param_set <- map_chr(1:nrow(param_comb), function(i) {
    paste0('l100', 'n', param_comb[i, 'n'], 'f', param_comb[i, 'f'], 'u', param_comb[i, 'u'])
})

####################################################
# check if all the simulated sequences are present #
####################################################

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

#####################################
# check if the Treeso result exists #
#####################################

for (param in param_set){
    for (ind in 1:run_size){
        treeso_file <- paste0(param, '/simresult_', param, '_', ind, '.txt')
        if (!file.exists(treeso_file)) {
            print(paste0('Missing: ', treeso_file))
        }
    }
}


#%% for the missing result, check the log files

missing_si 
