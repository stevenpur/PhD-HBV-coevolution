#%% import libraries
library(tidyverse)

#%% load data
sim_dir <- '~/hbv_covar3/analysis/sim_seq/' 
list.files(sim_dir, pattern = 'simresult.*')

for (i in 1:100){
  simresult <- read_tsv(paste0(sim_dir, 'simresult.', i, '.txt'))
  print(simresult %>% head())
}
