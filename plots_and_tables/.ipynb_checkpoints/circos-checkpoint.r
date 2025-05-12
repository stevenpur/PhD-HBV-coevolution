library(tidyverse)
library(circlize)


color_scheme <- list(
    # Major category: X
    X = "#A9A9A9",  # Dark Grey for the category X

    # Major category: P
    terminal_protein = "#FFA07A",      # Light Salmon
    spacer = "#FF8C00",                # Dark Orange
    reverse_transcriptase = "#FFD700", # Gold
    RNAse_h = "#FFC107",               # Orange Red

    # Major category: S
    pre_s1 = "#87CEEB",                # Sky Blue
    pre_s2 = "#4682B4",                # Steel Blue
    surface_antigen = "#1E90FF",       # Dodger Blue

    # Major category: C
    pre_core = "#98FB98",              # Pale Green
    core = "#2E8B57",                  # Sea Green

    # Repeated Colors for Major Categories
    P = "#FF8C00",                     # Dark Orange for category label
    S = "#1E90FF",                     # Dodger Blue for category label
    C = "#2E8B57"                      # Sea Green for category label
)


setwd("~/hbv_covar3")
genotypes <- c("A", "B", "C", "D")
genes <- c("C", "P", "S", "X")
potential_ovp_threshold <- 0

#%% load data
assoc_files <- paste0("./analysis/hbv", genotypes, "_lasso_nofuse_pair_extraInfo_mspval_hmp_domain.tsv")
names(assoc_files) <- genotypes

message("loading association files...")
assocs <- map(assoc_files, function(x) {
    read.table(x, header = T, stringsAsFactors = F)
})
names(assocs) <- genotypes

domain_pos_file <- "./analysis/protein_domain_pos.rds"
domain_pos <- readRDS(domain_pos_file)

#%% set up plotting df
set_plotting_df <- function(genotype){

    domain_data <- list()
    cur_y <- 0.3

    for (gene in genes) {
        cur_y <- cur_y + 0.5
        # we want to access the name of the domains, so we will use genotype A as an example
        for (domain in names(domain_pos[[genotype]][[gene]])) {
            start <- domain_pos[[genotype]][[gene]][[domain]][1]
            end <- domain_pos[[genotype]][[gene]][[domain]][2]
            if (end < start) {
                end <- end + 3664 + 20
            }
            domain_data[[domain]] <- c(domain, cur_y, cur_y + 0.5, start, end, color_scheme[[domain]])
            print(domain)
            print(domain_data[[domain]])
        }
    }
    domain_data <- data.frame(do.call(rbind, domain_data), stringsAsFactors = F)
    colnames(domain_data) <- c("name", "y_bottom", "y_top", "start", "end", "color")

    domain_data$y_bottom <- as.numeric(domain_data$y_bottom)
    domain_data$y_top <- as.numeric(domain_data$y_top)
    domain_data$start <- as.numeric(domain_data$start)
    domain_data$end <- as.numeric(domain_data$end)

    return(domain_data)
}

#%% plot
for (genotype in genotypes) {
    domain_data <- set_plotting_df(genotype)
    assoc <- assocs[[genotype]]
    assoc_sig <- assoc[which(assoc$adj_pval < 0.05 & assoc$potential_same_loci == 0 &
                            (assoc$bs1 > 0.9 | assoc$bs2 > 0.9)), ]

    # set up the basic circos plot
    options(repr.plot.width=10, repr.plot.height=10)
    pdf(paste0("~/hbv_covar3/plots/hbv", genotype, "_circos_plot.pdf"))
    genome_df <- data.frame(chr = '0', start = 1, end = 3664)
    circos.par(list('track.margin'=c(0.01,0.1)))
    circos.genomicInitialize(genome_df, tickLabelsStartFromZero = T, )
    circos.genomicRect(domain_data[,c("start", "end")], 'test',
                        ytop = domain_data$y_top,
                        ybottom = domain_data$y_bottom,
                        col = domain_data$color,)
    link_pos1 <- data.frame(chr = '0', 
                            start = assoc_sig$wg_start_pos1, 
                            end = assoc_sig$wg_start_pos1)
    link_pos2 <- data.frame(chr = '0', 
                            start = assoc_sig$wg_start_pos2, 
                            end = assoc_sig$wg_start_pos2)
    circos.genomicLink(region1 = link_pos1, region2 = link_pos2)
    dev.off()
}

