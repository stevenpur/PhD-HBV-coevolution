# simulate the tree
library(ape)

# get the population size
args <- commandArgs(trailingOnly = TRUE)
pop_size <- as.numeric(args[1])

# simulate the tree using coalescent model
tree <- rcoal(pop_size)

outfile <- paste0("/well/bag/clme1992/hbv_covar3/analysis/sim_seq/simseq_N", pop_size, ".tree")
write.tree(tree, outfile)