
## How the simulation was run

Supposing each site in the genome is biallelic, we define a rate matrix Q for a single site:

|     | A   | B   |
| --- | --- | --- |
| A   | -$u$  | $u$   |
| B   | $u$   | -$u$  | 

where $u$ is the mutation rate. Following this rate matrix, we can define a rate matrix $Q_{pair}$ for a pair of sites that are evolving independently:

|     | AA  | AB  | BA  | BB  |
| --- | --- | --- | --- | --- |
| AA  |  -2$u$   | $u$ | $u$ | 0   |
| AB  | $u$ |  -2$u$   | 0   | $u$ |
| BA  | $u$ | 0   |   -2$u$   | $u$ |
| BB  | 0   | $u$ | $u$ |    -2$u$  |

This matrix assumes that the probability of two sites changing simultaneously is 0. Now, if we consider a pair of sites that are coevolving. For a pair of sites (X, Y) that take the allele pair AA:

Site X evolves independently from Site Y, freely changing from A to B. Site Y is dependent on Site X, whenever Site X, changes from A to B, Site Y has an incentive to also change from A to B, and a disincentive to change from B to A.
   
   In such a scenario, we would get the following rate matrix $Q_{cov1}$
   
|     | AA  | AB  | BA  | BB  |
| --- | --- | --- | --- | --- |
| AA  |  -2$u$   | $u$ | $u$ | 0   |
| AB  | $u$ |  -2$u$   | 0   | $u$ |
| BA  | $u$ | 0   |   -$(1+C)u$   | $Cu$ |
| BB  | 0   | $u$ | $u/C$ |    -$(1+1/C)u$  |

Similarly, we can also consider 3 sites together (X, Y, Z) where X and Y have a coevolutionary relationship; Y and Z have a coevolutionary relationship. 

X -> Y -> Z

X is allowed to mutate freely, independent from Y and Z; Y is dependent on X but independent from Y; Z is dependent on Y, but independent from X. The matrix for these three sites is (skipping the diagonal element for simplicity):


|     | AAA   | ABA   | BAA       | BBA         | AAB | ABB  | BAB   | BBB  |
| --- | ----- | ----- | --------- | ----------- | --- | ---- | ----- | ---- |
| AAA | - | $u$   | $u$       | 0           | $u$ | 0    | 0     | 0    |
| ABA | $u$   | - | 0         | $u$         | 0   | $Cu$ | 0     | 0    |
| BAA | $u$   | 0     | - | $Cu$        | 0   | 0    | $u$   | 0    |
| BBA | 0     | $u$   | $u/C$     | -| 0   | 0    | 0     | $Cu$ |
| AAB | $u$   | 0     | 0         | 0           |  -   | $u$  | $u$   | 0    |
| ABB | 0     | $u/C$ | 0         | 0           | $u$ |  -    | 0     | $u$  |
| BAB | 0     | 0     | $u$       | 0           | $u$ | 0    |    -   | $Cu$ |
| BBB | 0     | 0     | 0         | $u/c$       | 0   | $u$  | $u/C$ |  -    | 

Now that we have the rate matrix, we can simulate the sequence accordingly:

1. simulate a tree

```R
tree <- rtree(pop_size, rooted = T)

# mimic the edge length distribution to the one we found in the actual HBV tree
hbv_tree <- read.tree("RAxML_bestTree.HBVA_withOutGroup_tree")
tree$edge.length <- sample(
						hbv_tree$edge.length, 
						length(tree$edge.length),
						replace = TRUE)
```

2. simulate the coevolving sites by using the simseq function. This function takes the parameter of 
	* a given tree
	* starting sequence at the root of the tree
	* The allele available (A and B for biallelic sites)
	* The rate matrix
	Then simulate the sequence by sampling from the $exp(Qt)$ matrix at each node of the tree. Where $t$ is the time passed (the branch length), and $Q$ is the rate matrix. 

	For example, for a given branch with length $t$, if at the start of the branch, a given site has an allele A. To determine what the state of the site is, we would need to sample the allele from the probability distribution: $X exp(Qt)$, where X is (1, 0) to represent allele A at the start of the branch.
	
	The full code of the function is as follow:
	
```R
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
```