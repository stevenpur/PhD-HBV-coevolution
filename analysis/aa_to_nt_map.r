###################################################################
# This code do two things:                                        #
# 1. Map AA msa location to NT msa location                       #
# 2. Identify the nucldeotides that correspond to each amino acid #
#    for the information of overlapping amino aicd in different   #
#    genes dureing downstream analysis                            #
###################################################################
library(tidyverse)
library(ape)
library(seqinr)
setwd("~/hbv_covar3")

# basic parameter
genotype <- "D"
genes <- c("C", "P", "S", "X")

# input file
gene_msa_pos_file <- "./QC/step4_translate_gene/gene_msa_pos.rds"
nt_msa_file <- "./QC/step4_translate_gene/HBV_ALLgeno_phyloFilter_ref_mafft.fasta"
aa_msa_files <- list()
for (gene in genes) {
    aa_msa_files[[gene]] <- paste0("./amino_acid/HBV_gene", gene, "_AAseqs_noShift_stopToDot_mafft.fasta")
}
# output file
wrt_file_aa_to_nt_pos <- paste0("./analysis/hbv", genotype, "_gene_aa_to_genome_pos_perSamp.rds")
wrt_file_aa_to_nt_id <- paste0("./analysis/hbv", genotype, "_gene_aa_to_genome_id_perSamp.rds")

# define function
skip_gap_get_AAcode <- function(nt_seq, pos) {
    cnt <- 0
    result <- vector(length = 3)
    nt_pos <- vector()
    while (cnt != 3) {
        if (pos > length(nt_seq)) {
            return(NULL)
        }
        if (nt_seq[pos] == "-") {
            pos <- pos + 1
            next
        } else {
            cnt <- cnt + 1
            result[cnt] <- nt_seq[pos]
            nt_pos <- c(nt_pos, pos)
            pos <- pos + 1
        }
    }
    return(list(nt_code = result, end_pos = pos - 1, nt_pos = nt_pos))
}

# read in NT msa and gene posistion
gene_msa_pos <- readRDS(gene_msa_pos_file)
nt_msa_seqs <- read.dna(nt_msa_file, format = "fasta", as.character = T, as.matrix = T)

gene_aa_to_genome_pos <- list() # the resulting map is stored here
gene_aa_to_nt_id <- list() # the result of nt ID for each aa is stored here

for (gene in genes) {
    message(gene)
    # wrangle AA msa seqs
    aa_msa_seqs <- seqinr::read.fasta(aa_msa_files[[gene]], seqtype = "AA")
    aa_msa_seqs <- do.call(rbind, aa_msa_seqs)
    aa_msa_genotypes <- map_chr(strsplit(rownames(aa_msa_seqs), "_"), 2)
    aa_msa_seqs <- aa_msa_seqs[which(aa_msa_genotypes == genotype), ]
    aa_msa_seqs[which(aa_msa_seqs == ".")] <- "*"
    aa_msa_seqs_cons <- seqinr::consensus(aa_msa_seqs)
    aa_msa_pos <- 1:(ncol(aa_msa_seqs))
    # wrangle NT msa seqs
    gene_nt_seqs <- nt_msa_seqs[, gene_msa_pos[[genotype]][[gene]]]
    gene_nt_pos <- gene_msa_pos[[genotype]][[gene]]
    gene_nt_genotypes <- map_chr(strsplit(rownames(gene_nt_seqs), "_"), 2)
    gene_nt_seqs <- gene_nt_seqs[which(gene_nt_genotypes == genotype), ]
    gene_nt_seqs_cons <- seqinr::consensus(gene_nt_seqs)
    # variables used to store informtion
    nt_ind <- matrix(nrow = nrow(aa_msa_seqs), ncol = ncol(aa_msa_seqs)) # used to record NT position for each MSA position
    rownames(nt_ind) <- rownames(aa_msa_seqs)
    cur_aa_to_nt_id <- matrix(nrow = nrow(aa_msa_seqs), ncol = ncol(aa_msa_seqs))
    rownames(cur_aa_to_nt_id) <- rownames(aa_msa_seqs)
    nt_seq <- vector()
    pos_start <- c(1)
    pos_end <- vector()
    for (i in seq_len(length(aa_msa_seqs_cons))) {
        message(paste0("site ", i))
        for (sample in rownames(aa_msa_seqs)) {
            last_nt_ind <- max(c(tail(na.exclude(nt_ind[sample, ]), 1), 1))
            get_nt_code_result <- skip_gap_get_AAcode(gene_nt_seqs[sample, ], last_nt_ind)
            sample_nt_code <- get_nt_code_result$nt_code
            if (aa_msa_seqs[sample, ][i] != "-") {
                nt_ind[sample, i] <- get_nt_code_result$end_pos + 1
                cur_aa_to_nt_id[sample, i] <- paste(gene_nt_pos[get_nt_code_result$nt_pos], collapse = '_')
                if (translate(sample_nt_code) != aa_msa_seqs[sample, i]) {
                    message("no gap case:")
                    message(sample)
                    message(paste0("gene ", gene, ":", i))
                    message(paste0(translate(sample_nt_code), " vs ", aa_msa_seqs[sample, ][i]))
                    stop("not equal 0")
                }
            }
        }
        nt_ind_freq <- table(nt_ind[, i])
        if (length(nt_ind_freq) > 1) {
            print(nt_ind_freq)
        }
        if (length(nt_ind_freq) == 0) {
            pos_start <- c(pos_start, NA)
            pos_end <- c(pos_end, NA)
            nt_seq <- c(nt_seq, NA)
        } else {
            matter_row <- which(nt_ind[, i] != "-")
            if (i != 1) {
                prev_matter_nt_ind_freq <- table(unlist(apply(nt_ind[matter_row, -i, drop = F], 1, function(x) {
                    tail(na.exclude(x), 1)
                })))
                pos_start <- c(pos_start, as.numeric(names(prev_matter_nt_ind_freq)[which.max(prev_matter_nt_ind_freq)]))
            }
            pos_end <- c(pos_end, as.numeric(names(nt_ind_freq)[which.max(nt_ind_freq)]) - 1)
            nt_seq <- c(nt_seq, paste(gene_nt_seqs_cons[tail(pos_start, 1):tail(pos_end, 1)], collapse = ""))
        }
    }
    wg_pos_start <- gene_nt_pos[pos_start]
    wg_pos_end <- gene_nt_pos[pos_end]
    data <- data.frame(
        aa_msa_pos = 1:ncol(aa_msa_seqs),
        codon_start = wg_pos_start,
        codon_end = wg_pos_end,
        nt_seq = nt_seq
    )
    gene_aa_to_genome_pos[[gene]] <- data
    gene_aa_to_nt_id[[gene]] <- cur_aa_to_nt_id
}

saveRDS(gene_aa_to_genome_pos, file = wrt_file_aa_to_nt_pos)
saveRDS(gene_aa_to_nt_id, file = wrt_file_aa_to_nt_id)