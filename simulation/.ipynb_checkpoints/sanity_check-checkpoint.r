#%% import libraries
library(tidyverse)

#%% check to see if all simulation result is available 
sim_dir <- '~/hbv_covar3/analysis/sim_seq/' 
simresult_files <- list.files(sim_dir, pattern = 'simresult.*')

param_set <- map_chr(sim_files, ~str_match(.x, 'simresult_(.*)\\_.*.txt')[,2])

param_set_unique <- unique(param_set)
param_set_cnt <- param_set_unique.map()
    vector('list', length(param_set_unique)) 
for (i in 1:length(param_set_unique)) {
    cur_param_set <- param_set_unique[i]
    pattern = paste0('_', cur_param_set, '_')
    cur_sim_files <- sim_fiels[str_detect(sim_files, pattern)]
    param_set_cnt[[i]] <- c(cur_param_set, as.character(length(cur_sim_files)))
}

#%% for the missing result, check the log files

missing_si 
