library(tidyverse)
library(seqinr)
library(foreach)
library(doParallel)
library(itertools)

setwd("~/hbv_covar3")
ncores <- 10
for (genotype in c("A", "B", "C", "D")) {
    message(genotype)
    # genotype <- "A"
    aa_to_nt_id_file <- paste0("./analysis/hbv", genotype, "_gene_aa_to_genome_id_perSamp.rds")
    aa_to_nt_id <- readRDS(aa_to_nt_id_file)
    allele_switch_file <- paste0("./analysis/hbv", genotype, "_allele_switch_cnt.rds")
    allele_switch <- readRDS(allele_switch_file)
    combination <- combn(names(allele_switch), 2)

    registerDoParallel(ncores)
    iterator_chunks <- isplitVector(1:ncol(combination), chunks = ncores)
    potential_same_loci <- foreach(chunk = iterator_chunks) %dopar% {
        result <- vector()
        for (i in chunk) {
            if (i %% 10000 == 0) {
                message(i)
            }
            site_id1 <- combination[1, i]
            site_id2 <- combination[2, i]
            gene1 <- gsub("gene", "", gsub("_site.*", "", site_id1))
            gene2 <- gsub("gene", "", gsub("_site.*", "", site_id2))
            if (gene1 == gene2) {
                result <- c(result, 0)
                next
            }
            aa_pos1 <- as.numeric(gsub(".*_site", "", site_id1))
            aa_pos2 <- as.numeric(gsub(".*_site", "", site_id2))
            id1 <- strsplit(aa_to_nt_id[[gene1]][, aa_pos1], "_")
            id2 <- strsplit(aa_to_nt_id[[gene2]][, aa_pos2], "_")
            samp_names <- intersect(names(id1), names(id2))
            id1 <- id1[samp_names]
            id2 <- id2[samp_names]
            na_samp_ind <- which(is.na(id1) & is.na(id2))
            if (length(na_samp_ind) != 0) {
                id1 <- id1[-na_samp_ind]
                id2 <- id2[-na_samp_ind]
            }
            id_overlap <- map_dbl(seq_len(length(id1)), function(x) {
                overlap_len <- length(na.exclude(intersect(id1[[x]], id2[[x]])))
                if (overlap_len > 0) {
                    return(1)
                } else {
                    return(0)
                }
            })
            result <- c(result, sum(id_overlap) / length(id1))
        }
        return(result)
    }
    potential_same_loci <- unlist(potential_same_loci)
    ovp_matrix <- matrix(ncol = length(allele_switch), nrow = length(allele_switch))
    colnames(ovp_matrix) <- names(allele_switch)
    rownames(ovp_matrix) <- names(allele_switch)
    for (i in 1:length(potential_same_loci)) {
        site_id1 <- combination[1, i]
        site_id2 <- combination[2, i]
        ovp_matrix[site_id1, site_id2] <- potential_same_loci[i]
        ovp_matrix[site_id2, site_id1] <- potential_same_loci[i]
    }
    saveRDS(ovp_matrix, paste0("~/hbv_covar3/analysis/hbv", genotype, "_overlap.rds"))
}