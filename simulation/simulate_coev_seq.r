# This script simulates sequences with coevolving sites given a tree.

# Required Libraries:
# - ape: for tree manipulation and simulation
# - phangorn: for phylogenetic analysis
# - Matrix: for matrix operations
# - tidyverse: for data manipulation

# Arguments:
# - --len: length of the simulated sequences
# - --tree: file path to the tree file
# - --coev_factor: coevolution factor for coevolving site pairs
# - --run_ind: run index for multiple simulations
# - --mu: mutation rate

# The script checks, searches, and parses the command line arguments.
# It writes the simulation result to the output file in FASTA format.

message('loading libraries...')
library(ape)
library(phangorn)
library(Matrix)
library(tidyverse)

message('loading user parameters')
# arguments: --len, --tree, --coev_factor, --run_ind, --mu
# check, search, and parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (!("--len" %in% args)) {
    stop("Missing argument: --len")
}
if (!("--tree" %in% args)) {
    stop("Missing argument: --tree")
}
if (!("--coev_factor" %in% args)) {
    stop("Missing argument: --coev_factor")
}
if (!("--run_ind" %in% args)) {
    stop("Missing argument: --run_ind")
}
if (!("--mu" %in% args)) {
    stop("Missing argument: --mu")
}

seq_len <- as.numeric(args[which(args == "--len") + 1])
tree_file <- as.character(args[which(args == "--tree") + 1])
coev_factor <- as.numeric(args[which(args == "--coev_factor") + 1])
run_ind <- as.numeric(args[which(args == "--run_ind") + 1])
u <- as.numeric(args[which(args == "--mu") + 1])

message('loading tree files...')
# get the population size from the tree file
tree <- read.tree(tree_file)
pop_size <- length(tree$tip.label)
# set the run id
run_id <- paste0("l", seq_len, "n", pop_size, "f", coev_factor, "u", u, "_", run_ind)

message("seq_len: ", seq_len, "\n",
        "tree_file: ", tree_file, "\n",
        "pop_size: ", pop_size, "\n",
        "coev_factor: ", coev_factor, "\n",
        "run_ind: ", run_ind, "\n",
        "u: ", u, "\n")

# set the output file
outfile <- paste0("~/hbv_covar3/analysis/sim_seq/test/simseq_", run_id, ".fasta")

# set bases of the sequences, assuming biallelic for all sites
bases <- c("x", "y")
coev_pair <- c("yy", "xx")

message("setting Q matrices...")
# set Q of independant sites
Q <- matrix(c(-1*u, u, u, -1*u), nrow = 2, ncol = 2)
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
siteAB_coev_pair <- c("yy", "xx")
siteBC_coev_pair <- c("yy", "xx")
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

message("Q matrix:")
print(Qco3)


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
    node_name <- paste0('iid', cur_node)
    result[[node_name]]  <- start_seq
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
            # has not reached tip yet, go down another level
            child_result <- simseq(tree, child_seq, levs, Q, child)
            # concat child result to the current result
            result <- c(result, child_result)
        } else {
            node_name <- tree$tip.label[child]
            result[[node_name]] <- child_seq
        }
    }
    return(result)
}

sim_seq_ind <- function(tree, rate) {
    # function implementation

    # Create a comprehensive state array for all nodes
    all_node_states <- rep(0, Nnode(tree) + length(tree$tip.label))

    # Calculate total branch length
    total_branch_length <- sum(tree$edge.length)

    # Step 1: Determine number of substitution events
    num_events <- rpois(1, lambda = rate * total_branch_length)

    # Step 2: Determine the location of substitution events
    event_locations <- runif(num_events, min = 0, max = total_branch_length)

    # Initialize the alleles at all nodes (0 = ancestral allele, 1 = derived allele)
    tree$node.label <- rep(0, Nnode(tree))

    # Function to propagate change to descendant nodes
    propagateChange <- function(tree, node, newState, all_node_states) {
        all_node_states[node] <- newState
        descendants <- Descendants(tree, node, type = "all")
        all_node_states[descendants] <- newState
        return(all_node_states)
    }

    # Simulate allele changes at each event
    current_length <- 0
    node_ind <- 0
    for (event in sort(event_locations)) {
        while (current_length < event) {
            node_ind <- node_ind + 1
            current_length <- current_length + tree$edge.length[node_ind]
            
        }
        node <- tree$edge[node_ind, 2]
        newState <- 1 - all_node_states[node]
        all_node_states <- propagateChange(tree, node, newState, all_node_states)
    }

    # Step 4: Retrieve the allele each node
    alleles <- all_node_states
    names(alleles) <- c(tree$tip.label, paste0("iid", length(tree$tip.label)+(1:Nnode(tree))))
    return(alleles)
}


# simulate the coevolving sites
message("simulating coevolving sites...")
start_time <- Sys.time()
sim_result_coev <- simseq(tree, c("xxx"), base_tri, Qco3)
end_time <- Sys.time()
message("time used: " , end_time - start_time)

# put the simulation result into a matrix
sim_msa <- do.call(rbind, strsplit(unlist(sim_result_coev), ""))


# simulate the independant sites
message("simulating independant independant sites...")
start_time <- Sys.time()
sim_ind_result <- list()
for (site in 1:seq_len) {
    sim_ind_result[[site]] <- ifelse(sim_seq_ind(tree, u) == 1, "x", "y")
}
sim_ind_result <- do.call(cbind, sim_ind_result)
end_time <- Sys.time()
message("time used: " , end_time - start_time)

# merge the independant sites to the coevolving sites by row names
sim_msa <- cbind(sim_msa, sim_ind_result[rownames(sim_msa), ])
# write the simulation result
write.dna(sim_msa, outfile, format = "fasta", colsep = "")
