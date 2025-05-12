library(tidyverse)
library(seqinr)
library(grid)
library(tidyverse)
library(patchwork)

setwd('~/hbv_covar3')

genotype_color <- list(
     A = '#F8766D',
     B = '#7CAE00',
     C = '#00BFC4', 
     D = '#C77CFF'
)

genes <- c("S", "P", "X", "C")
genotypes <- c('A', 'B', 'C', 'D')

assoc_files <- paste0("./analysis/hbv", genotypes, "_lasso_nofuse_pair_extraInfo_mspval_hmp_domain.tsv")
names(assoc_files) <- genotypes

assocs <- map(assoc_files, function(x) {
    read.table(x, header = T, stringsAsFactors = F)
})
names(assocs) <- genotypes


get_assoc_dist_info <- function(assoc, dist_3d_tbl){
    assoc_dist <- vector()
    assoc_aa_dist <- vector()
    assoc_genotype <- vector()
    assoc_name <- vector()
    no_assoc_dist <- dist_3d_tbl
    for (genotype in genotypes){
        assoc <- assocs[[genotype]]
        assoc_sig <- assoc[which(assoc$adj_pval < 0.05 &
                                 (assoc$bs1 > 0.9 | assoc$bs2 > 0.9)), ]
        for(i in 1:nrow(assoc_sig)) {
            siteA <- assoc_sig$siteA[i]
            siteB <- assoc_sig$siteB[i]
            if (all(c(siteA, siteB) %in% rownames(dist_3d_tbl))){
                cur_aa_dist <- abs(assoc_sig$gene_aa_pos1[i] - assoc_sig$gene_aa_pos2[i])
                assoc_aa_dist <- c(assoc_aa_dist, cur_aa_dist)
                assoc_genotype <- c(assoc_genotype, genotype)
                assoc_name <- c(assoc_name, assoc_sig$site_pair_id[i])
                assoc_dist <- c(assoc_dist, dist_3d_tbl[siteA, siteB])
                no_assoc_dist[siteA, siteB] <- NA
            }
        }
    }
    no_assoc_dist <- na.exclude(no_assoc_dist[lower.tri(no_assoc_dist, diag=FALSE)])
    return (list(assoc_dist = assoc_dist,
                 assoc_aa_dist = assoc_aa_dist,
                 assoc_genotype = assoc_genotype, 
                 assoc_name = assoc_name,
                 no_assoc_dist = no_assoc_dist))
}

######## for core protein #########
core_3d_tb <- read.table("./analysis/geneC_structure_distance_matrix.tsv")
core_dist_info <- get_assoc_dist_info(assocs, core_3d_tb)
t.test(core_dist_info$assoc_dist, core_dist_info$no_assoc_dist, var.equal = F)

# plotting
data <- as.data.frame(rbind(cbind(core_dist_info$assoc_dist, "covarying"), cbind(core_dist_info$no_assoc_dist, "non-covarying")))
colnames(data) <- c("distance", "type")
data$distance <- as.numeric(data$distance)
core_box_plot <- ggplot(data, aes(x = type, y = distance)) + 
    geom_boxplot() +
    theme(axis.title=element_blank(),
          text = element_text(size=20))

ggsave("./plots/core_3D_dist.png", core_box_plot)


####### for polymerase protein #########
p_3d_tb <- read.table("./analysis/geneP_structure_distance_matrix.tsv")
p_dist_info <- get_assoc_dist_info(assocs, p_3d_tb)
t.test(p_dist_info$assoc_dist, p_dist_info$no_assoc_dist, var.equal = F)

# plotting
data <- as.data.frame(rbind(cbind(p_dist_info$assoc_dist, "covarying"), cbind(p_dist_info$no_assoc_dist, "non-covarying")))
colnames(data) <- c("distance", "type")
data$distance <- as.numeric(data$distance)
p_box_plot <- ggplot(data, aes(x = type, y = distance)) + 
    geom_boxplot() +
    theme(axis.title=element_blank(),
          text = element_text(size=20))

ggsave("./plots/polymerase_3D_dist.png", p_box_plot)

###########################
# stratify by aa distance #
###########################
get_dist_byaa <- function(assoc_aa_dists, dist_3d_tb, gene, exclude_pairs = c()){
    dist3d_byaa <- list()
    sites_with_3d <- as.numeric(gsub(".*site", "", names(dist_3d_tb)))
    for (aa_dist in assoc_aa_dists){
        start_aa <- 1
        end_aa <- 1 + aa_dist
        cur_3d_dist <- vector()
        while(end_aa <= max(sites_with_3d)){
            siteA <- paste0("gene", gene, "_site", start_aa)
            siteB <- paste0("gene", gene, "_site", end_aa)
            if (paste(sort(c(siteA, siteB)), collapse = "_") %in% exclude_pairs){
                start_aa <- start_aa + 1
                end_aa <- end_aa + 1
                next
            }
            if (all(c(siteA, siteB) %in% names(dist_3d_tb))){
                cur_3d_dist <- c(dist_3d_tb[siteA, siteB], cur_3d_dist)
            }
            start_aa <- start_aa + 1
            end_aa <- end_aa + 1
        }
        dist3d_byaa[[paste0("dist", aa_dist)]] <- cur_3d_dist
    }
    return(dist3d_byaa)
}

plot_strat_violin <- function(data){
    p <- ggplot(data, aes(x = aa_dist, y = dist3d)) +
        geom_violin() + 
        geom_point(data = data[which(data$type == "dist"),], 
                   position = position_jitter(seed = 1, h = 0),
                   alpha = 0.2,
                   size = 0.4) +
        geom_point(aes(colour = I(color), shape = assoc_name), 
                    data = data[which(data$type == "assoc"),], 
                    position = position_jitter(seed = 1, h = 0), 
                    size = 3) +
        theme(legend.position="none",
              axis.title=element_blank(),
              text = element_text(size=20))
    return(p)
}


c_assoc_aa_dists <- core_dist_info$assoc_aa_dist
c_dist3d_byaa <- get_dist_byaa(c_assoc_aa_dists, core_3d_tb, "C")

p_assoc_aa_dists <- p_dist_info$assoc_aa_dist
p_dist3d_byaa <- get_dist_byaa(p_dist_info$assoc_aa_dist, p_3d_tb, "P")

p_dist3d_byaa_extra <- get_dist_byaa(c(152), p_3d_tb, "P")

plot_dist3d_byaa <- function(dist3d_byaa, dist_info, plt_nrows = 2){
    plots <- list()
    for (i in 1:length(dist3d_byaa)){
        aa_dist <- as.numeric(gsub("dist", "", names(dist3d_byaa)[i]))
        cur_dist_ind <- which(dist_info$assoc_aa_dist == aa_dist)
        cur_assoc_dist <- dist_info$assoc_dist[cur_dist_ind]
        cur_genotype <- dist_info$assoc_genotype[cur_dist_ind]
        cur_name <- dist_info$assoc_name[cur_dist_ind]
        data <- rbind(cbind(aa_dist, cur_assoc_dist, "assoc", cur_genotype, cur_name),
                          cbind(aa_dist, dist3d_byaa[[i]], "dist", NA, NA))

        data <- as.data.frame(data)
        colnames(data) <- c("aa_dist", "dist3d", "type", "genotype", "assoc_name")
        data$dist3d <- as.numeric(data$dist3d)
        data$color <- map_chr(data$genotype, function(x) {
            ifelse(is.na(x), "black", genotype_color[[x]])
        })
        data$aa_dist <- paste0(data$aa_dist, " aa apart")
        p <- plot_strat_violin(data)
        plots[[i]] <- p
    }
    myplot <- patchwork::wrap_plots(plots, nrow=plt_nrows)
    return(myplot)
}

core_strat_plot <- plot_dist3d_byaa(c_dist3d_byaa, core_dist_info, 2)
p_strat_plot <- plot_dist3d_byaa(p_dist3d_byaa, p_dist_info, 1)

extra_info <- list(
    'assoc_aa_dist' = c(152),
    'assoc_dist'= c(p_3d_tb['geneP_site469', 'geneP_site621']),
    'assoc_genotype' = 'B',
    'assoc_name' = 'geneP_site469_geneP_site621'
)

p_strat_extra_plot <- plot_dist3d_byaa(p_dist3d_byaa_extra, extra_info, 1)




ggsave("./plots/core_3D_dist_byaa.png", p, width = 14, height = 7)

combine_plot <- (plot_spacer() + core_box_plot + p_box_plot + plot_spacer() + plot_layout(widths = c(1, 2, 2, 1))) / 
    plot_spacer() /
    (core_strat_plot) /
    plot_spacer() /
    (plot_spacer() + p_strat_plot + plot_spacer() + p_strat_extra_plot + plot_spacer() + plot_layout(widths = c(0.5, 3, 0.5, 1, 0.5))) +
    plot_layout(heights = c(1.5, 0.25, 2, 0.25, 1)) + theme(plot.margin=margin(100,0,100,0))

ggsave('./plots/dist_combined.png', combine_plot, width = 14, height = 14)
