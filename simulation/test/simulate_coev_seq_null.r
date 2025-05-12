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
if (!("--run_ind" %in% args)) {
    stop("Missing argument: --run_ind")
}
if (!("--mu" %in% args)) {
    stop("Missing argument: --mu")
}

seq_len <- as.numeric(args[which(args == "--len") + 1])
tree_file <- as.character(args[which(args == "--tree") + 1])
run_ind <- as.numeric(args[which(args == "--run_ind") + 1])
u <- as.numeric(args[which(args == "--mu") + 1])

message('loading tree files...')
# get the population size from the tree file
tree <- read.tree(tree_file)
pop_size <- length(tree$tip.label)
# set the run id
param_set <- paste0("l", seq_len, "n", pop_size, "u", u)
run_id <- paste0(param_set, "_", run_ind)

message("seq_len: ", seq_len, "\n",
        "tree_file: ", tree_file, "\n",
        "pop_size: ", pop_size, "\n",
        "run_ind: ", run_ind, "\n",
        "u: ", u, "\n")


# set the output file
outfile <- paste0("~/hbv_covar3/analysis/sim_seq/", param_set, "/simseq_", run_id, "rescaled.fasta")

# set bases of the sequences, assuming biallelic for all sites
bases <- c("x", "y")
coev_pair <- c("yy", "xx")

message("setting Q matrices...")
# set Q of independant sites
Q <- matrix(c(-1*u, u, u, -1*u), nrow = 2, ncol = 2)
colnames(Q) <- bases
rownames(Q) <- bases

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


result <- sim_ind_result[tree$tip.label, ]

# write the simulation result
write.dna(result, outfile, format = "fasta", colsep = "")
