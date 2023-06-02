library(ape)

setwd("~/hbv_covar3")
genotypes <- c("A", "B", "C", "D")

for (genotype in genotypes) {
    # output file
    wrt_file <- paste0("./analysis/phylo_build/RAxML_bestTree.HBV", genotype, "_rooted")

    # input files
    tree_file <- paste0("./analysis/phylo_build/RAxML_bestTree.HBV", genotype, "_withOutGroup_tree")
    tree <- read.tree(tree_file)

    # root the tree and remove outgroup
    out_group <- grep("ref", tree$tip.label)
    tree <- root(tree, outgroup = out_group)
    tree <- drop.tip(tree, out_group)

    write.tree(tree, wrt_file)
}