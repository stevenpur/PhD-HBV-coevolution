library(tidyverse)
library(ggplot2)
setwd("~/hbv_covar3/analysis/sim_seq/")
# make a plot where the y axis is whether the coevolution is identified or not, the x-axis is the number of co-mutation in the tree?
 
# first, we need to identify the number of mutations in the tree 

#%% define functions
get_lnfu_comb <- function(files) {
    files <- gsub(".*/", "", files)
    files_strat <- map_chr(strsplit(files, "_"), 2)
    files_strat <- strsplit(files_strat, '[lnfu]')
strat_tb <- do.call(rbind, files_strat)[, -1]
strat_tb <- as.data.frame(apply(strat_tb, 2, as.numeric))
    colnames(strat_tb) <- c('l', 'n', 'f', 'u')
    l_uniq <- unique(strat_tb$l)
    n_uniq <- unique(strat_tb$n)
    f_uniq <- unique(strat_tb$f)
    u_uniq <- unique(strat_tb$u)
    lnfu_comb <- expand.grid(l_uniq, n_uniq, f_uniq, u_uniq)
    colnames(lnfu_comb) <- c('l', 'n', 'f', 'u')
    lnfu_comb$param_set <- paste0('l', lnfu_comb$l, 'n', lnfu_comb$n, 'f', lnfu_comb$f, 'u', lnfu_comb$u)
    return(lnfu_comb)
}

mycolors <- function(x) {
   colors <- colorRampPalette(c("darkblue", "green", "yellow"))(8)
}

#Function to create the contour plot
create_contour_plot <- function(data, z, mybreaks, mycolors, title, position) {
   p <- ggplot(data, aes(x=log10(coev_factor), y=log10(mutation_rate), z=z)) +
      geom_contour_filled(breaks = mybreaks) +
      scale_fill_manual(name = title, values = mycolors(10), drop=FALSE) +
      scale_x_continuous(name="Log of Coevolution Factor (f)", 
                         breaks=log10(data$coev_factor), 
                         labels=data$coev_factor) +
      scale_y_continuous(name="Log of Mutation Rate (u)", 
                         breaks=log10(data$mutation_rate), 
                         labels=data$mutation_rate) +
      theme(legend.position = "none") +
      ggtitle(title)
   # if position is not left, then set the text to ""
   if (position != "left") {
      p <- p + theme(axis.title.y = element_text(color = "transparent"))
   }
   # if position is not middle, then hide the plot title (by setting it to empty string)
   if (position != "middle") {
      p <- p + theme(plot.title = element_text(color = "transparent"))
   }
   # if position is not right then hide the legend
   if (position == "legend") {
      p <- p + theme(legend.position = "right")
      p_legend <- ggplot_gtable(ggplot_build(p))
      legend <- p_legend$grobs[[which(p_legend$layout$name == "guide-box")]]
      return(legend)
   }
   return(p)
}

get_tp_fp <- function(file, coev_pair, indir_coev_pair, l){
    if (file.info(file)$size < 2) {
        tb <- NA 
    } else {
        tb <- read.csv(file, header = F, sep = " ")
        colnames(tb) <- c("siteA", "siteB")
        tb$pair_id <- map_chr(1:nrow(tb), function(x) {
            paste(sort(c(tb[x, 1], tb[x, 2])), collapse = ":")
        })
    }
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
    return (result)
}

get_site_switch <- function(log_file){
    log <- readRDS(log_file)
    site1_switch <- log$site_switch$site1
    site2_switch <- log$site_switch$site2
    site3_switch <- log$site_switch$site3
    site12_co_switch <- log$site_switch$site1 & log$site_switch$site2
    site23_co_switch <- log$site_switch$site2 & log$site_switch$site3
    site13_co_switch <- log$site_switch$site1 & log$site_switch$site3
    result <- c(sum(site1_switch, na.rm = TRUE),
                sum(site2_switch, na.rm = TRUE),
                sum(site3_switch, na.rm = TRUE),
                sum(site12_co_switch, na.rm = TRUE),
                sum(site23_co_switch, na.rm = TRUE), 
                sum(site13_co_switch, na.rm = TRUE))
    names(result) <- c("site1", "site2", "site3", "site1_2", "site2_3", "site1_3")
    return(result)
}



#%% plot the result of the simulation
coev_pair <- c("site1:site2", "site2:site3")
indir_coev_pair <- c("site1:site3")

# files <- list.files('test/', pattern = "log_l.*.rds", full.names = TRUE)
pattern <- 'l100'
folders <- list.dirs(full.names = FALSE, recursive = FALSE)
folders <- folders[str_detect(folders, pattern)]

lnfu_comb <- get_lnfu_comb(files)
u_uniq <- unique(lnfu_comb$u)
f_uniq <- unique(lnfu_comb$f)
n_uniq <- unique(lnfu_comb$n)

l <- 100
n <- 100

comb_fu <- expand.grid(f_uniq, u_uniq)
colnames(comb_fu) <- c('f', 'u')
data <- list()

for (i in 1:nrow(comb_fu)){
    f <- comb_fu$f[i]
    u <- comb_fu$u[i]
    param_set <- paste0('l100n', n, 'f', f, 'u', u)
    result_files <- list.files(param_set, pattern = paste0('simresult_', param_set, '_.*.txt'), full.names = TRUE)
    cur_data <- map(result_files, function(file){
        sim_ind <- gsub('.txt', '', strsplit(file, "_")[[1]][3])
        log_file <- paste0('./', param_set, '/log_', param_set, '_', sim_ind, '.rds')
        tp_fp <- get_tp_fp(file, coev_pair, indir_coev_pair, l)
        site_switch <- get_site_switch(log_file)
        result <- c(param_set, sim_ind, tp_fp, site_switch)
        names(result)[1] <- "param_set"
        names(result)[2] <- "sim_ind"
        return(result)
    })
    data[[length(data) + 1]] <- do.call(rbind, cur_data)
}

data_df <- as.data.frame(do.call(rbind, data))
for (i in 2:ncol(data_df)){
    data_df[, i] <- as.numeric(data_df[, i])
}

# i have the data, but they are distributed too densely
# I need to bin the data
data_df <- data_df %>%
    mutate(bin = ntile(site1_2, 5))
plt_df <- data_df
plt_df$bin <- as.factor(plt_df$bin)

p <- ggplot(plt_df, aes(x=bin, y=tp)) +
    geom_jitter(alpha = 0.1, width = 0.2, height = 0.1)

ggsave("~/hbv_covar3/plots/comutation_vs_tp.png", p, width=8, height=6, units='in')

plt_df <- as.data.frame(do.call(rbind, result_summary))
colnames(plt_df) <- c('coev_factor', 'mutation_rate', 'co_mut_cnt')

p <- ggplot(plt_df, aes(x=log10(coev_factor), y=log10(mutation_rate), z=co_mut_cnt)) +
      geom_contour_filled()

# save the plot
ggsave("~/hbv_covar3/plots/coev_factor_vs_mutation_rate.png", p, width=8, height=6, units="in")


