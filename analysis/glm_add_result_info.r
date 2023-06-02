library(tidyverse)
library(seqinr)
library(foreach)
library(doParallel)
library(itertools)

# basic parameter
args <- commandArgs(trailingOnly = TRUE)

setwd("~/hbv_covar3/analysis/")
genotype <- "B"
genes <- c("C", "P", "S", "X")
ncores <- 12

# get the association result
HBV.assoc <- list()
file_name <- paste0("./hbv", genotype, "_glm_edge_test_lrt.tsv")
HBV.assoc <- read.table(file_name, sep = "\t", stringsAsFactors = F, header = T)
# read in aaMSA to ntMSA position map
gene_aa_to_genome_pos <- readRDS(paste0("./hbv", genotype, "_gene_aa_to_genome_pos_perSamp.rds"))
# read in nt's id for each amino acid
aa_to_nt_id_file <- paste0("./hbv", genotype, "_gene_aa_to_genome_id_perSamp.rds")
aa_to_nt_id <- readRDS(aa_to_nt_id_file)

# add genome position info column
wg_data <- apply(HBV.assoc, 1, function(x) {
    gene1 <- gsub("gene", "", gsub("_site.*", "", x[["siteA"]]))
    gene2 <- gsub("gene", "", gsub("_site.*", "", x[["siteB"]]))
    data1 <- gene_aa_to_genome_pos[[gene1]]
    data2 <- gene_aa_to_genome_pos[[gene2]]
    aa_pos1 <- as.numeric(gsub(".*site", "", x[["siteA"]]))
    aa_pos2 <- as.numeric(gsub(".*site", "", x[["siteB"]]))
    row1 <- which(data1$aa_msa_pos == aa_pos1)
    row2 <- which(data2$aa_msa_pos == aa_pos2)
    wg_start1 <- data1$codon_start[row1]
    wg_start2 <- data2$codon_start[row2]
    wg_end1 <- data1$codon_end[row1]
    wg_end2 <- data2$codon_end[row2]
    return(c(gene1, gene2, aa_pos1, aa_pos2, wg_start1, wg_start2, wg_end1, wg_end2))
})

wg_data <- t(wg_data)
HBV.assoc$gene1 <- wg_data[, 1]
HBV.assoc$gene2 <- wg_data[, 2]
HBV.assoc$gene_aa_pos1 <- as.numeric(wg_data[, 3])
HBV.assoc$gene_aa_pos2 <- as.numeric(wg_data[, 4])
HBV.assoc$wg_start_pos1 <- as.numeric(wg_data[, 5])
HBV.assoc$wg_start_pos2 <- as.numeric(wg_data[, 6])
HBV.assoc$wg_end_pos1 <- as.numeric(wg_data[, 7])
HBV.assoc$wg_end_pos2 <- as.numeric(wg_data[, 8])

# add potential same location info column
HBV.assoc.sort <- HBV.assoc[order(HBV.assoc$site_pair_id), ]
registerDoParallel(ncores)
iterator_chunks <- isplitVector(seq_len(nrow(HBV.assoc.sort)), chunks = ncores)
potential_same_loci <- foreach(chunk = iterator_chunks) %dopar% {
    result <- vector()
    cur_site_pair_id <- ""
    for (i in chunk) {
        if (i %% 10000 == 0) {
            message(i)
        }
        if (HBV.assoc.sort$site_pair_id[i] == cur_site_pair_id) {
            result <- c(result, result[length(result)])
            next
        } else {
            cur_site_pair_id <- HBV.assoc.sort$site_pair_id[i]
        }
        gene1 <- HBV.assoc.sort$gene1[i]
        gene2 <- HBV.assoc.sort$gene2[i]
        if (gene1 == gene2) {
            result <- c(result, 0)
            next
        }
        aa_pos1 <- HBV.assoc.sort$gene_aa_pos1[i]
        aa_pos2 <- HBV.assoc.sort$gene_aa_pos2[i]
        id1 <- strsplit(aa_to_nt_id[[gene1]][, aa_pos1], "_")
        id2 <- strsplit(aa_to_nt_id[[gene2]][, aa_pos2], "_")
        samp_names <- intersect(names(id1), names(id2))
        id1 <- id1[samp_names]
        id2 <- id2[samp_names]
        id_overlap <- map_dbl(seq_len(length(samp_names)), function(x) {
            overlap_len <- length(na.exclude(intersect(id1[[x]], id2[[x]])))
            if (overlap_len > 0) {
                return(1)
            } else {
                return(0)
            }
        })
        result <- c(result, sum(id_overlap) / length(samp_names))
    }
    return(result)
}
HBV.assoc.sort$potential_same_loci <- unlist(potential_same_loci)


# add fdr adjusted pval
HBV.assoc.sort$pval_fdr <- p.adjust(HBV.assoc.sort$pval, method = "fdr")

# write result
wrt_file <- paste0("./hbv", genotype, "_glm_edge_test_extraInfo.tsv")
write.table(HBV.assoc.sort, wrt_file, sep = "\t", quote = F, col.names = T, row.names = F)