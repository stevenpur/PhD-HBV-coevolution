# simulate sequences with coevolving sites given a tree #

library(ape)
library(phangorn)
library(Matrix)
library(tidyverse)
library(ggtree)

args <- commandArgs(trailingOnly = TRUE)
run_id <- args[1]
outfile <- paste0("/users/bag/hlq763/hbv_covar3/analysis/sim_seq/simseq_", run_id, ".txt")
outfile_tree <- paste0("/users/bag/hlq763/hbv_covar3/analysis/sim_seq/simseq_", run_id, ".tree")

# set simulation parameter, assuming biallelic sites (allele x and y)
seq_len <- 100
pop_size <- 1000
bases <- c("x", "y")
coev_pair <- c("yy")
coev_factor <- 100

# set Q of independant sites
Q <- matrix(c(-1, 1, 1, -1), nrow = 2, ncol = 2)
colnames(Q) <- bases
rownames(Q) <- bases

# set Q of coevolving site pairs
base_pairs <- expand.grid(bases, bases)
base_pairs <- paste0(base_pairs[, 1], base_pairs[, 2])
Qco <- matrix(nrow = 4, ncol = 4)
colnames(Qco) <- base_pairs
rownames(Qco) <- base_pairs

# siteA is independantly evolving, siteB is following siteA
for (i in 1:nrow(Qco)) {
    for (j in 1:ncol(Qco)) {
        # skip if there is no change on either sites
        if (i == j) {
            next
        }
        siteA_change <- c(strsplit(base_pairs[i], "")[[1]][1], strsplit(base_pairs[j], "")[[1]][1])
        siteB_change <- c(strsplit(base_pairs[i], "")[[1]][2], strsplit(base_pairs[j], "")[[1]][2])

        # treat changes in both sites simultaneously as impossible to happen
        if (siteA_change[1] != siteA_change[2] && siteB_change[1] != siteB_change[2]) {
            Qco[i, j] <- 0
            next
        }

        # treat one site changes as the changes in the independant model
        # first, get the changing base
        base_change <- NA
        if (siteA_change[1] != siteA_change[2]) {
            base_change <- siteA_change
            # if the change is at siteA, and not at siteB, then this element is not influence by coev pair because siteA evolve idependantly
            Qco[i, j] <- Q[base_change[1], base_change[2]]
            next
        } else {
            base_change <- siteB_change
        }
        # get the rate of this chanage
        Qco[i, j] <- Q[base_change[1], base_change[2]]

        # now, consider if coevolution should be consider
        if (base_pairs[i] %in% coev_pair) {
            Qco[i, j] <- Qco[i, j] / coev_factor
        } else if (base_pairs[j] %in% coev_pair) {
            Qco[i, j] <- Qco[i, j] * coev_factor
        }
    }
}

for (i in 1:nrow(Qco)) {
    Qco[i, i] <- -1 * sum(Qco[i, ], na.rm = T)
}

# three site coevolution model
Qco3 <- matrix(nrow = 8, ncol = 8)
base_tri <- expand.grid(bases, bases, bases)
base_tri <- paste0(base_tri[, 1], base_tri[, 2], base_tri[, 3])
colnames(Qco3) <- base_tri
rownames(Qco3) <- base_tri
siteAB_coev_pair <- "yy"
siteBC_coev_pair <- "yy"
for (i in 1:nrow(Qco3)) {
    for (j in 1:ncol(Qco3)) {
        # skip if there is no change
        if (i == j) {
            next
        }
        siteA_change <- c(strsplit(base_tri[i], "")[[1]][1], strsplit(base_tri[j], "")[[1]][1])
        siteB_change <- c(strsplit(base_tri[i], "")[[1]][2], strsplit(base_tri[j], "")[[1]][2])
        siteC_change <- c(strsplit(base_tri[i], "")[[1]][3], strsplit(base_tri[j], "")[[1]][3])

        siteA_is_change <- as.numeric(siteA_change[1] != siteA_change[2])
        siteB_is_change <- as.numeric(siteB_change[1] != siteB_change[2])
        siteC_is_change <- as.numeric(siteC_change[1] != siteC_change[2])
        # treat changes in >2 sites simultaneously as impossible to happen
        if (siteA_is_change + siteB_is_change + siteC_is_change >= 2) {
            Qco3[i, j] <- 0
            next
        }
        # consider siteA change, this site is changing in independantly from other sites
        if (siteA_is_change == 1) {
            Qco3[i, j] <- Q[siteA_change[1], siteA_change[2]]
            next
        }
        # consider siteB change, this site is dependant on the second site
        if (siteB_is_change == 1) {
            # first set the change as if siteB is changing independantly
            Qco3[i, j] <- Q[siteB_change[1], siteB_change[2]]
            # get the bases that siteA and B are changing from
            baseAB <- paste0(siteA_change[1], siteB_change[1])
            if (baseAB %in% siteAB_coev_pair) {
                Qco3[i, j] <- Qco3[i, j] / coev_factor
            }
            # get the bases that siteA and B are changing into
            baseAB <- paste0(siteA_change[2], siteB_change[2])
            if (baseAB %in% siteAB_coev_pair) {
                Qco3[i, j] <- Qco3[i, j] * coev_factor
            }
            next
        }
        # consider siteC change, this site is dependant on the second site
        if (siteC_is_change == 1) {
            # first set the change as if siteC is changing independantly
            Qco3[i, j] <- Q[siteC_change[1], siteC_change[2]]
            # get the bases that siteB and C are changing from
            baseBC <- paste0(siteB_change[1], siteC_change[1])
            if (baseBC %in% siteBC_coev_pair) {
                Qco3[i, j] <- Qco3[i, j] / coev_factor
            }
            # get the bases that siteB and C are changing into
            baseBC <- paste0(siteB_change[2], siteC_change[2])
            if (baseBC %in% siteBC_coev_pair) {
                Qco3[i, j] <- Qco3[i, j] * coev_factor
            }
            next
        }
    }
}


for (i in 1:nrow(Qco3)) {
    Qco3[i, i] <- -1 * sum(Qco3[i, ], na.rm = T)
}




# function for sequence simulation

simseq <- function(tree, start_seq, levs, Q, cur_node = length(tree$tip.label) + 1, result = list(), head_run = T) {
    # get the starting state for each base in the sequence
    seq_states <- map(start_seq, function(cur_base) {
        as.numeric(cur_base == levs)
    })
    seq_states <- do.call(cbind, seq_states) # each column is the state of a base
    children <- Children(tree, cur_node)
    for (child in children) {
        # get branch length
        branch_len <- tree$edge.length[which(tree$edge[, 2] == child)]
        child_seq <- map_chr(1:ncol(seq_states), function(i) {
            # get the probability of the state after time t (the parent branch length)
            cur_state <- seq_states[, i]
            state_dist <- cur_state %*% expm(branch_len * Q)
            # get the base of the child based on the distribution
            cur_base <- sample(levs, 1, prob = as.vector(state_dist))
            return(cur_base)
        })
        if (child > length(tree$tip.label)) {
            # has not reach tip yet, go down another level
            result <- simseq(tree, child_seq, levs, Q, child, result, F)
        } else {
            node_name <- tree$tip.label[child]
            result[[node_name]] <- child_seq
        }
    }
    if (head_run) {
        # get the result into table format before output
        result_names <- names(result)
        result <- do.call(rbind, result)
        rownames(result) <- result_names
        result[match(tree$tip.label, rownames(result)), ]
    }
    return(result)
}


# generate a tree for simulation
tree <- rtree(pop_size, rooted = T)
# pool edge length from real-world HBV tree
hbv_tree <- read.tree("~/hbv_covar3/analysis/phylo_build/RAxML_bestTree.HBVA_withOutGroup_tree")
tree$edge.length <- sample(hbv_tree$edge.length, length(tree$edge.length), replace = TRUE)

# simulate the coevolving sites
# sim_result_coev <- simseq(tree, c("xx"), base_pairs, Qco)
sim_result_coev <- simseq(tree, c("xxx"), base_tri, Qco3)
# collect the simulatoin result of each tip of the tree
sim_msa <- do.call(rbind, strsplit(sim_result_coev[, 1], ""))

# simulate the independant sites
sim_result_ind <- simseq(tree, sample(bases, seq_len, replace = T), bases, Q)
sim_msa <- cbind(sim_msa, sim_result_ind)

# write the simulation result
write.table(sim_msa, outfile, row.names = TRUE, col.names = FALSE, quote = FALSE, sep = " ")
write.tree(tree, outfile_tree)
