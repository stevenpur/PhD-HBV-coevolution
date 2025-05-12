# This script simulates sequences with coevolving sites given a tree.

# Required Libraries:
# - ape: for tree manipulation and simulation
# - phangorn: for phylogenetic analysis
# - Matrix: for matrix operations
# - tidyverse: for data manipulation
# - optparse: for command line argument parsing
# - logging: for proper logging

# Arguments:
# - --len: length of the simulated sequences
# - --tree: file path to the tree file
# - --coev_factor: coevolution factor for coevolving site pairs
# - --run_ind: run index for multiple simulations
# - --mu: mutation rate
# - --covar_desc_per: percentage of descendants for co-mutation events
# - --output_dir: output directory for simulation results

# The script checks, searches, and parses the command line arguments.
# It writes the simulation result to the output file in FASTA format.

# Load required libraries
message('Loading libraries...')
suppressPackageStartupMessages({
    library(ape)
    library(phangorn)
    library(Matrix)
    library(tidyverse)
    library(optparse)
    library(logging)
})

# Configuration
config <- list(
    bases = c("x", "y"),
    coev_pairs = c("yy", "xx")
)

# Setup logging
basicConfig(level = 'INFO')
loginfo("Starting sequence simulation")

# Parse command line arguments
option_list <- list(
    make_option(c("--len"), type = "integer", help = "Length of the simulated sequences"),
    make_option(c("--tree"), type = "character", help = "File path to the tree file"),
    make_option(c("--coev_factor"), type = "double", help = "Coevolution factor for coevolving site pairs"),
    make_option(c("--run_ind"), type = "integer", help = "Run index for multiple simulations"),
    make_option(c("--mu"), type = "double", help = "Mutation rate"),
    make_option(c("--covar_desc_per"), type = "double", default = NULL, help = "Optional: Percentage of descendants for co-mutation events"),
    make_option(c("--output_dir"), type = "character", default = "./", help = "Output directory for simulation results")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate arguments
validate_args <- function(opt) {
    if (is.null(opt$len) || opt$len <= 0) {
        stop("Sequence length must be a positive integer")
    }
    if (is.null(opt$tree) || !file.exists(opt$tree)) {
        stop("Tree file must exist and be readable")
    }
    if (is.null(opt$coev_factor) || opt$coev_factor <= 0) {
        stop("Coevolution factor must be positive")
    }
    if (is.null(opt$run_ind) || opt$run_ind < 0) {
        stop("Run index must be non-negative")
    }
    if (is.null(opt$mu) || opt$mu <= 0) {
        stop("Mutation rate must be positive")
    }
    if (!is.null(opt$covar_desc_per) && (opt$covar_desc_per < 0 || opt$covar_desc_per > 1)) {
        stop("Co-mutation percentage must be between 0 and 1")
    }
    if (!dir.exists(opt$output_dir)) {
        stop("Output directory must exist")
    }
}

tryCatch({
    validate_args(opt)
}, error = function(e) {
    logerror("Argument validation failed: %s", e$message)
    stop(e)
})

# Set parameters
seq_len <- opt$len
tree_file <- opt$tree
coev_factor <- opt$coev_factor
run_ind <- opt$run_ind
u <- opt$mu

# Load and validate tree
loginfo("Loading tree file: %s", tree_file)
tree <- tryCatch({
    read.tree(tree_file)
}, error = function(e) {
    logerror("Failed to read tree file: %s", e$message)
    stop(e)
})

pop_size <- length(tree$tip.label)
param_set <- paste0("l", seq_len, "n", pop_size, "f", coev_factor, "u", u)
run_id <- paste0(param_set, "_", run_ind)

# Log parameters
loginfo("Simulation parameters:")
loginfo("Sequence length: %d", seq_len)
loginfo("Population size: %d", pop_size)
loginfo("Coevolution factor: %f", coev_factor)
loginfo("Run index: %d", run_ind)
loginfo("Mutation rate: %f", u)

outfile <- file.path(opt$output_dir, paste0("simseq_", run_id, "rescaled.fasta"))

# Matrix setup functions
setup_rate_matrix <- function(u) {
    # Set up the rate matrix for independent sites
    # This matrix represents the transition rates between states x and y
    # -u represents the rate of leaving a state
    # u represents the rate of entering a state
    loginfo("Setting up rate matrix with mutation rate: %f", u)
    Q <- matrix(c(-1*u, u, u, -1*u), nrow = 2, ncol = 2)
    colnames(Q) <- config$bases
    rownames(Q) <- config$bases
    return(Q)
}

setup_coev_matrix <- function(Q, coev_factor) {
    # Set up the rate matrix for coevolving site pairs
    # This matrix represents the transition rates between all possible pairs of states (xx, xy, yx, yy)
    # The coevolution factor affects the rates when sites are in coevolving states
    loginfo("Setting up coevolution matrix with factor: %f", coev_factor)
    base_pairs <- expand.grid(config$bases, config$bases)
    base_pairs <- paste0(base_pairs[, 1], base_pairs[, 2])
    Qco <- matrix(nrow = 4, ncol = 4)
    colnames(Qco) <- base_pairs
    rownames(Qco) <- base_pairs

    for (i in 1:nrow(Qco)) {
        for (j in 1:ncol(Qco)) {
            if (i == j) next
            
            # Get the state changes for both sites
            siteA_change <- c(strsplit(base_pairs[i], "")[[1]][1], strsplit(base_pairs[j], "")[[1]][1])
            siteB_change <- c(strsplit(base_pairs[i], "")[[1]][2], strsplit(base_pairs[j], "")[[1]][2])

            # Treat changes in both sites simultaneously as impossible
            if (siteA_change[1] != siteA_change[2] && siteB_change[1] != siteB_change[2]) {
                Qco[i, j] <- 0
                next
            }

            # Get the changing base and its rate from the independent model
            base_change <- if (siteA_change[1] != siteA_change[2]) siteA_change else siteB_change
            Qco[i, j] <- Q[base_change[1], base_change[2]]

            # Apply coevolution factor:
            # - If starting from a coevolving state, reduce the rate of leaving it
            # - If ending in a coevolving state, increase the rate of entering it
            if (base_pairs[i] %in% config$coev_pairs) {
                Qco[i, j] <- Qco[i, j] / coev_factor
            } else if (base_pairs[j] %in% config$coev_pairs) {
                Qco[i, j] <- Qco[i, j] * coev_factor
            }
        }
    }

    # Set diagonal elements to negative sum of row (conservation of probability)
    for (i in 1:nrow(Qco)) {
        Qco[i, i] <- -1 * sum(Qco[i, ], na.rm = TRUE)
    }
    
    return(Qco)
}

setup_three_site_matrix <- function(Q, coev_factor, u) {
    # Set up the rate matrix for three coevolving sites
    # This matrix represents the transition rates between all possible triplets of states (xxx, xxy, etc.)
    # Sites A and B are coevolving, as are sites B and C
    loginfo("Setting up three-site coevolution matrix")
    base_tri <- expand.grid(config$bases, config$bases, config$bases)
    base_tri <- paste0(base_tri[, 1], base_tri[, 2], base_tri[, 3])
    Qco3 <- matrix(nrow = 8, ncol = 8)
    colnames(Qco3) <- base_tri
    rownames(Qco3) <- base_tri

    for (i in 1:nrow(Qco3)) {
        for (j in 1:ncol(Qco3)) {
            if (i == j) next
            
            # Get the state changes for all three sites
            siteA_change <- c(strsplit(base_tri[i], "")[[1]][1], strsplit(base_tri[j], "")[[1]][1])
            siteB_change <- c(strsplit(base_tri[i], "")[[1]][2], strsplit(base_tri[j], "")[[1]][2])
            siteC_change <- c(strsplit(base_tri[i], "")[[1]][3], strsplit(base_tri[j], "")[[1]][3])

            # Check which sites are changing
            changes <- c(
                siteA_change[1] != siteA_change[2],
                siteB_change[1] != siteB_change[2],
                siteC_change[1] != siteC_change[2]
            )
            
            # Treat changes in more than one site simultaneously as impossible
            if (sum(changes) >= 2) {
                Qco3[i, j] <- 0
                next
            }

            # Site A changes independently
            if (changes[1]) {
                Qco3[i, j] <- Q[siteA_change[1], siteA_change[2]]
                next
            }

            # Site B changes, affected by coevolution with site A
            if (changes[2]) {
                # Start with independent rate
                Qco3[i, j] <- Q[siteB_change[1], siteB_change[2]]
                # Apply coevolution factor based on states of A and B
                baseAB <- paste0(siteA_change[1], siteB_change[1])
                if (baseAB %in% config$coev_pairs) {
                    Qco3[i, j] <- Qco3[i, j] / coev_factor
                }
                baseAB <- paste0(siteA_change[2], siteB_change[2])
                if (baseAB %in% config$coev_pairs) {
                    Qco3[i, j] <- Qco3[i, j] * coev_factor
                }
                next
            }

            # Site C changes, affected by coevolution with site B
            if (changes[3]) {
                # Start with independent rate
                Qco3[i, j] <- Q[siteC_change[1], siteC_change[2]]
                # Apply coevolution factor based on states of B and C
                baseBC <- paste0(siteB_change[1], siteC_change[1])
                if (baseBC %in% config$coev_pairs) {
                    Qco3[i, j] <- Qco3[i, j] / coev_factor
                }
                baseBC <- paste0(siteB_change[2], siteC_change[2])
                if (baseBC %in% config$coev_pairs) {
                    Qco3[i, j] <- Qco3[i, j] * coev_factor
                }
            }
        }
        Qco3 <- normalise_Q(Qco3, u)
    }

    # Set diagonal elements to negative sum of row (conservation of probability)
    for (i in 1:nrow(Qco3)) {
        Qco3[i, i] <- -1 * sum(Qco3[i, ], na.rm = TRUE)
    }
    
    return(Qco3)
}

# Simulation functions
simulate_sequence <- function(tree, start_seq, levs, Q, cur_node = length(tree$tip.label) + 1, result = list()) {
    # Simulate sequence evolution along a phylogenetic tree
    # For each site:
    # 1. Convert the current state to a probability vector
    # 2. Calculate the probability distribution after the branch length
    # 3. Sample the new state from this distribution
    tryCatch({
        # Convert sequence to state probability vectors
        seq_states <- map(start_seq, function(cur_base) {
            as.numeric(cur_base == levs)
        })
        seq_states <- do.call(cbind, seq_states)
        
        node_name <- paste0('iid', cur_node)
        result[[node_name]] <- start_seq
        
        # Process each child node
        children <- Children(tree, cur_node)
        for (child in children) {
            # Get branch length and calculate new states
            branch_len <- tree$edge.length[which(tree$edge[, 2] == child)]
            child_seq <- map_chr(1:ncol(seq_states), function(i) {
                cur_state <- seq_states[, i]
                # Calculate probability distribution after branch length
                state_dist <- cur_state %*% expm(branch_len * Q)
                # Sample new state
                sample(levs, 1, prob = as.vector(state_dist))
            })
            
            # Recursively process child nodes
            if (child > length(tree$tip.label)) {
                child_result <- simulate_sequence(tree, child_seq, levs, Q, child)
                result <- c(result, child_result)
            } else {
                node_name <- tree$tip.label[child]
                result[[node_name]] <- child_seq
            }
        }
        return(result)
    }, error = function(e) {
        logerror("Error in sequence simulation: %s", e$message)
        stop(e)
    })
}

simulate_independent_sites <- function(tree, rate) {
    # Simulate independent site evolution using a continuous-time Markov process
    # 1. Calculate total branch length and number of events
    # 2. Place events randomly along the tree
    # 3. Propagate changes to descendant nodes
    tryCatch({
        all_node_states <- rep(0, Nnode(tree) + length(tree$tip.label))
        total_branch_length <- sum(tree$edge.length)
        num_events <- rpois(1, lambda = rate * total_branch_length)
        event_locations <- runif(num_events, min = 0, max = total_branch_length)
        
        tree$node.label <- rep(0, Nnode(tree))
        
        # Helper function to propagate state changes to all descendants
        propagateChange <- function(tree, node, newState, all_node_states) {
            all_node_states[node] <- newState
            descendants <- Descendants(tree, node, type = "all")
            all_node_states[descendants] <- newState
            return(all_node_states)
        }
        
        # Process each event
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
        
        alleles <- all_node_states
        names(alleles) <- c(tree$tip.label, paste0("iid", length(tree$tip.label)+(1:Nnode(tree))))
        return(alleles)
    }, error = function(e) {
        logerror("Error in independent site simulation: %s", e$message)
        stop(e)
    })
}

# Add after the validate_args function
mutate_xy <- function(nucleotide) {
    if (nucleotide == 'x') {
        return('y')
    } else if (nucleotide == 'y') {
        return('x')
    } else {
        stop(paste('unexpected nucleotide:', nucleotide))
    }
}

get_stationary_distribution <- function(Q) {
    # Calculate the stationary distribution of a rate matrix
    # This is the eigenvector corresponding to the eigenvalue 0
    # The stationary distribution is the long-term equilibrium distribution of the Markov process
    # It satisfies the equation Q * pi = 0
    ev <- eigen(t(Q))
    # Find the eigenvalue closest to zero
    idx0 <- which.min(abs(ev$values))
    # Extract the corresponding (complex) eigenvector, take real part
    v <- Re(ev$vectors[, idx0])
    # Normalize to sum to 1
    pi <- v / sum(v)
    # Clip any tiny negative numerical artefacts
    pi[pi < 0] <- 0
    return(pi / sum(pi))
}

normalise_Q <- function(Q, u) {
    # Normalize the rate matrix to have the average rate of u
    Q <- Q / sum(Q) * u
    return(Q)
}

scale_Q_to_u <- function(Q, u) {
    # 1. stationary distribution of Q
    pi <- get_stationary_distribution(Q)
    # 2. current average (per‐state) rate
    R_current <- sum(pi * (-diag(Q)))
    # 3. scale factor to make average = u
    factor <- u / R_current
    # 4. return rescaled matrix
    Q * factor
}

# Add before the main simulation section
# Select a node for co-mutation
desc_n <- map_dbl(1:(length(tree$tip.label) + Ntip(tree)), function(x) {
    length(unlist(Descendants(tree, x, type="tips")))
})

covar_desc_per <- opt$covar_desc_per
loop_cnt <- 0
covar_nodes <- NULL
covar_desc_size <- round(pop_size * covar_desc_per)

while (length(covar_nodes) == 0) {
    covar_nodes <- which(desc_n == covar_desc_size-loop_cnt)
    loop_cnt <- loop_cnt + 1
}
if (length(covar_nodes) > 1) {
    covar_node <- sample(covar_nodes, 1)
} else {
    covar_node <- covar_nodes
}

descendants <- Descendants(tree, covar_node, type = "all")
descendants_labels <- map_chr(descendants, function(x) {
    if (x <= length(tree$tip.label)){
        return(tree$tip.label[x])
    } else {
        return(paste0('iid', x))
    }
})

# Main simulation
loginfo("Setting up rate matrices...")
Q <- setup_rate_matrix(u)
Qco <- setup_coev_matrix(Q, coev_factor)
Qco3 <- setup_three_site_matrix(Q, coev_factor, u)

# Get base_tri from Qco3 matrix
base_tri <- colnames(Qco3)

loginfo("Starting coevolving sites simulation...")
start_time <- Sys.time()
sim_result_coev <- simulate_sequence(tree, c("xxx"), base_tri, Qco3)
end_time <- Sys.time()
loginfo("Coevolving sites simulation completed in %f seconds", as.numeric(end_time - start_time))

loginfo("Starting independent sites simulation...")
start_time <- Sys.time()
sim_ind_result <- list()
for (site in 1:seq_len) {
    sim_ind_result[[site]] <- ifelse(simulate_independent_sites(tree, u) == 1, "x", "y")
}
sim_ind_result <- do.call(cbind, sim_ind_result)
end_time <- Sys.time()
loginfo("Independent sites simulation completed in %f seconds", as.numeric(end_time - start_time))

# Combine results
sim_msa <- do.call(rbind, strsplit(unlist(sim_result_coev), ""))
sim_msa <- cbind(sim_msa, sim_ind_result[rownames(sim_msa), ])

# Apply co-mutation if specified
if (!is.null(opt$covar_desc_per)) {
    loginfo("Applying co-mutation to sites 4 and 5 for %.2f%% of descendants", opt$covar_desc_per * 100)
    
    # Select a node for co-mutation
    desc_n <- map_dbl(1:(length(tree$tip.label) + Ntip(tree)), function(x) {
        length(unlist(Descendants(tree, x, type="tips")))
    })

    covar_desc_per <- opt$covar_desc_per
    loop_cnt <- 0
    covar_nodes <- NULL
    covar_desc_size <- round(pop_size * covar_desc_per)

    while (length(covar_nodes) == 0) {
        covar_nodes <- which(desc_n == covar_desc_size-loop_cnt)
        loop_cnt <- loop_cnt + 1
    }
    if (length(covar_nodes) > 1) {
        covar_node <- sample(covar_nodes, 1)
    } else {
        covar_node <- covar_nodes
    }

    descendants <- Descendants(tree, covar_node, type = "all")
    descendants_labels <- map_chr(descendants, function(x) {
        if (x <= length(tree$tip.label)){
            return(tree$tip.label[x])
        } else {
            return(paste0('iid', x))
        }
    })

    sim_msa_comut <- sim_msa
    comut_sites <- c(4, 5)
    for (site in comut_sites) {
        for (descendant in descendants_labels){
            sim_msa_comut[descendant, site] <- mutate_xy(sim_msa[descendant, site])
        }
    }
    final_msa <- sim_msa_comut
} else {
    loginfo("No co-mutation specified, using standard coevolution simulation")
    final_msa <- sim_msa
}

# Write results
loginfo("Writing results to %s", outfile)
tryCatch({
    write.dna(final_msa, outfile, format = "fasta", colsep = "")
    loginfo("Simulation completed successfully")
}, error = function(e) {
    logerror("Failed to write output file: %s", e$message)
    stop(e)
})
