library(tidyverse)

setwd("~/hbv_covar3")
genotypes <- c("A", "B", "C", "D")
genes <- c("C", "P", "S", "X")

#%% load the data

assoc_files <- paste0("./analysis/hbv", genotypes, "_lasso_nofuse_pair_extraInfo_mspval_hmp_domain.tsv")
names(assoc_files) <- genotypes
message("loading association files...")
assocs <- map(assoc_files, function(x) {
    read.table(x, header = T, stringsAsFactors = F)
})

#%% set the domains information
ovp_domains <- list(
    c("X", "pre_core"),
    c("X", "RNAse_h"),
    c("core", "terminal_protein"),
    c("pre_s1", "spacer"),
    c("pre_s2", "spacer"),
    c("pre_s2", "reverse_transcriptase"),
    c("surfice_antigen", "reverse_transcriptase")
)

domains <- c("terminal_protein", "spacer", 'reverse_transcriptase', 'RNAse_h',
             'pre_s1', 'pre_s2', 'surface_antigen', 
             'pre_core', 'core',
             'X')

domains <- unique(unlist(ovp_domains))
ovp_domains <- map_chr(ovp_domains, function(x) paste(sort(c(x[1], x[2])), collapse = ":"))                  



#%% plot the data
domain_rename <- list(
    terminal_protein="Terminal Protein",
    spacer = "Spacer",
    reverse_transcriptase = "Reverse Transcriptase",
    RNAse_h = "RNAse H",
    pre_s1 = "Pre-S1",
    pre_s2 = "Pre-S2",
    surface_antigen = "Surface Antigen",
    pre_core = "Pre-Core",
    core = "Core",
    X = "X"
)

plots <- list()
for (genotype in genotypes){
    assoc <- assocs[[genotype]]
    assoc <- assoc[which(assoc$potential_same_loci == 0), ]
    assoc_sig_cnt <- sum(assoc$adj_pval < 0.05 & 
                         (assoc$bs1 > 0.9 | assoc$bs2 > 0.9))
    assoc_nonsig_cnt <- nrow(assoc) - assoc_sig_cnt
    sig_p <- assoc_sig_cnt/nrow(assoc)
    heat_data <- data.frame()
    for (i in 1:length(domains)) {
        for (j in 1:i) {
            cur_domain_pair <- paste(sort(c(domains[i], domains[j])), collapse = ":")
            cur_assoc <- filter(assoc, domain_pair == cur_domain_pair)
            cur_sig_cnt <- sum(cur_assoc$adj_pval < 0.05 &
                              (cur_assoc$bs1 > 0.9 | cur_assoc$bs2 > 0.9))
            cur_nonsig_cnt <- nrow(cur_assoc) - cur_sig_cnt
            fisher_t <- matrix(c(cur_sig_cnt, assoc_sig_cnt-cur_sig_cnt, 
                              cur_nonsig_cnt, assoc_nonsig_cnt-cur_nonsig_cnt), 
                               nrow = 2, ncol = 2)
            fisher_test <- fisher.test(fisher_t, alternative = "greater")
            OR = fisher_test$estimate
            pval = fisher_test$p.value
            
            is_domain_ovp <- ifelse(cur_domain_pair %in% ovp_domains, 1, 0)
            heat_data <- rbind(heat_data, c(domains[i], domains[j], pval, OR, is_domain_ovp, cur_sig_cnt))
        }
    }
    colnames(heat_data) <- c("domain1", "domain2", "p_val", "OR", "is_ovp", "sig_cnt")
    heat_data$adj_pval <- p.adjust(heat_data$p_val, method = 'bonferroni')
    for (i in 1:nrow(heat_data)){
        heat_data <- rbind(heat_data, c(heat_data[i, 2],
                                        heat_data[i, 1],
                                        heat_data[i, 3],
                                        heat_data[i, 4],
                                        heat_data[i, 5],
                                        heat_data[i, 6],
                                        heat_data[i, 7]))
    }
    heat_data$OR <- as.numeric(heat_data$OR)
    heat_data$p_val <- as.numeric(heat_data$p_val)
    heat_data$sig_cnt <- as.numeric(heat_data$sig_cnt)
    heat_data$adj_pval <- as.numeric(heat_data$adj_pval)
    heat_data$is_ovp <- as.numeric(heat_data$is_ovp)
    heat_data$is_sig <- ifelse(heat_data$adj_pval < 0.05, 1, NA)
    heat_data$domain1 <- factor(map_chr(heat_data$domain1, function(x) domain_rename[[x]]),
                                levels = map_chr(domains, function(x) domain_rename[[x]]))
    heat_data$domain2 <- factor(map_chr(heat_data$domain2, function(x) domain_rename[[x]]),
                                levels = map_chr(domains, function(x) domain_rename[[x]]))
    heat_data$r <- heat_data$sig_cnt/max(heat_data$sig_cnt) * 5
    p <- ggplot(heat_data, aes(x=domain1, y=domain2)) + 
        geom_tile(fill = "white") +
        geom_point(aes(size = sig_cnt, colour = OR)) +
        scale_size_continuous(range = c(0, 10)) +
        scale_color_gradient2(low = "navajowhite1", high = "firebrick", mid = "navajowhite1") +
        geom_point(aes(x=domain1, y=domain2, size = is_sig), shape = "*") +
        #geom_tile(data = heat_data[which(heat_data$is_ovp == 1), ], 
        #          aes(x=domain1, y=domain2),
        #          fill = "transparent", 
        #          colour = "black",
        #          size = 1) +
        theme(axis.title = element_blank(), 
              axis.text.x = element_text(angle = 30, hjust=1)) +
        ggtitle(paste0("Genotype ", genotype)) +
        labs(size = "No. covaring residues")
    plots[[genotype]] <- p
}

options(repr.plot.width = 13, repr.plot.height = 10)
myplot <- patchwork::wrap_plots(plots, nrow=2)
print(myplot)
ggsave(filename="./plots/domains_heatmap.pdf", plot=myplot, width=13, height=10)

