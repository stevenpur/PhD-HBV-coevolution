library(tidyverse)
library(igraph)
library(ggraph)

setwd("~/hbv_covar3")
genotypes <- c("A", "B", "C", "D")
potential_ovp_threshold <- 0
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


#%% load files
assoc_files <- paste0("./analysis/hbv", genotypes, "_lasso_nofuse_pair_extraInfo_mspval_hmp_domain.tsv")
names(assoc_files) <- genotypes

assocs <- map(assoc_files, function(x) {
    read.table(x, header = T, stringsAsFactors = F)
})
names(assocs) <- genotypes

#%% plot graph
nets <- list()
for (genotype in genotypes){
    assoc <- assocs[[genotype]]
    # get the significant associations, their domains, and their strength (bootstrap value)
    assoc_sig <- assoc[which(assoc$adj_pval < 0.05 & 
                             assoc$potential_same_loci == 0 &
                             (assoc$bs1 > 0.9 | assoc$bs2 > 0.9)), ]

    assoc_sig <- assoc_sig %>% mutate(siteA_domain = paste0(siteA, "_d:", domain1),
                              siteB_domain = paste0(siteB, "_d:", domain2))

    weights <- map_dbl(1:nrow(assoc_sig), function(x) {
        max(assoc_sig$bs1[x], assoc_sig$bs2[x])
    })

    # create a graph from the significant
    graph_data <- assoc_sig[, c("siteA_domain", "siteB_domain")]
    graph_data$weight <- weights
    colnames(graph_data) <- c("nodeA", "nodeB", "weight")
    nodes <- unique(c(graph_data$nodeA, graph_data$nodeB))
    node_doms <- gsub(".*:", "", nodes)
    node_doms_colors <- map_chr(node_doms, function(cur_gene) color_scheme[[cur_gene]])
    net <- graph_from_data_frame(graph_data, vertices = nodes, directed = F)
    nets[[genotype]] <- net

    # plot the graph
    p <- ggraph(net, layout = 'fr') + 
        geom_edge_link(edge_color = "black", alpha = 0.7, width = 2) +  # Plot edges with transparency
        geom_node_point(aes(color = node_doms_colors), size = 8) +  # Plot nodes
        scale_color_identity() +  # Use identity scale for node colors
        theme_void()  # Remove axes and background
    ggsave(paste0("~/hbv_covar3/plots/hbv", genotype, "_assoc_graph2.pdf"), p, width = 10, height = 10)
}
