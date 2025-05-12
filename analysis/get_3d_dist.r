library(tidyverse)
library(Rpdb)
library(ape)
library(seqinr)
library(grid)
library(tidyverse)


pdb_extract_3d_tbl <- function(pdb, chain = '', elename = ''){
    result <- pdb$atoms
    if (chain != ''){
        result <- result[which(result$chainid == chain), ]
    }
    if (elename != ''){
        result <- result[which(result$elename == elename), ]
    }
}

add_siteid_to_3d_tbl <- function(tbl, msa, dom_start, gene) {
    msa_cons <- seqinr::consensus(msa)
    msa_pos <- 1:length(msa_cons)
    msa_pos_nogap <- msa_pos[-which(msa_cons == "-")]
    dom_start_nogap <- which(msa_pos_nogap == dom_start)
    tbl$site_id <- map_chr(tbl$resid, function(resid) {
        nogap_pos <- dom_start_nogap - 1 + resid
        gap_pos <- msa_pos_nogap[nogap_pos]
        return(paste0("gene", gene, "_site", gap_pos))
    })
    return(tbl)
}

eucl_dist <- function(v1, v2){
    return(sqrt(sum((v1-v2)^2)))
}

cal_residue_dist <- function(pos_tb, single_chain = FALSE) {
    if (single_chain) {
        rownames(pos_tb) <- pos_tb$site_id
        pos_tb <- pos_tb[, c("x1", "x2", "x3")]
        rownames(pos_tb) <- pos_tb$site_id
        dist_tbl <- as.matrix(dist(pos_tb))
        return(dist_tbl)
    }
    residue_num <- max(pos_tb$resid)
    dist_tbl <- matrix(nrow = residue_num, ncol = residue_num)
    rownames(dist_tbl) <- unique(pos_tb$site_id)
    rownames(dist_tbl) <- unique(pos_tb$site_id)
    for (i in 1:residue_num){
        site_id1 <- rownames(dist_tbl)[i]
        row1s <- which(pos_tb$site_id == site_id1)
        for (j in 1:residue_num){
            site_id2 <- rownames(dist_tbl)[j]
            row2s <- which(pos_tb$site_id == site_id2)
            row_comb <- expand.grid(row1s, row2s)
            cur_dist <- map_dbl(1:nrow(row_comb), function(x) {
                row1 <- row_comb[x, 1]
                row2 <- row_comb[x, 2]
                eucl_dist(c(pos_tb$x1[row1], pos_tb$x2[row1], pos_tb$x3[row1]),
                          c(pos_tb$x1[row2], pos_tb$x2[row2], pos_tb$x3[row2]))
            })
            dist_tbl[i,j]<- min(cur_dist)
        }
    }
    return(dist_tbl)
}

##### set parameters #####
genes <- c("S", "P", "X", "C")
genotypes <- c('A', 'B', 'C', 'D')
setwd("~/hbv_covar3")

######### read data ###############

gene_msa_files <- map(genes, function(x) {
    paste0("./amino_acid/HBV_gene", x, "_AAseqs_noShift_stopToDot_mafft.fasta")
})
names(gene_msa_files) <- genes
gene_msa <- map(genes, function(x) {
    do.call(rbind, read.fasta(gene_msa_files[[x]]))
})
names(gene_msa) <- genes

c_pdb <- read.pdb("./raw_data/7abl.pdb")
p_pdb <- read.pdb("./raw_data/pone.0106324.s001.pdb")

domain_pos_file <- "./analysis/protein_domain_pos_aa.rds"
domain_pos <- readRDS(domain_pos_file)$A

####################################


########## for core protein ################

core_dom_start <- domain_pos$C$core[1]
c_pdb_tb <- pdb_extract_3d_tbl(c_pdb, elename = 'CA')
c_pdb_tb <- add_siteid_to_3d_tbl(c_pdb_tb, gene_msa$C, core_dom_start, "C")
dm <- cal_residue_dist(c_pdb_tb)

write.table(dm, "./analysis/geneC_structure_distance_matrix.tsv")

############################################


######### for polymerase protein ###########

rt_dom_start <- domain_pos$P$reverse_transcriptase[1]
p_pdb_tb <- pdb_extract_3d_tbl(p_pdb, elename = 'CA')
p_pdb_tb <- add_siteid_to_3d_tbl(p_pdb_tb, gene_msa$P, rt_dom_start, "P")
dm <- cal_residue_dist(p_pdb_tb, single_chain = TRUE)

