bases <- c("x", "y")
coev_pair <- c("yy", "xx")

message("setting Q matrices...")
# set Q of independant sites
Q <- matrix(c('-1u', 'u', 'u', '-1*u'), nrow = 2, ncol = 2)
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
            Qco[i, j] <- '0'
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
            Qco[i, j] <- paste0(Qco[i, j], ' / f'
        } else if (base_pairs[j] %in% coev_pair) {
            Qco[i, j] <- paste0(Qco[i, j], 'f') 
        }
    }
}

for (i in 1:nrow(Qco)) {
    Qco[i, i] <- -1 * sum(Qco[i, ], na.rm = T)
}

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
                Qco3[i, j] <- paste0(Qco3[i, j],  '/f') 
            }
            # get the bases that siteA and B are changing into
            baseAB <- paste0(siteA_change[2], siteB_change[2])
            if (baseAB %in% siteAB_coev_pair) {
                Qco3[i, j] <- paste0(Qco3[i, j], 'f')
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
                Qco3[i, j] <- paste0(Qco3[i, j] , 'f')
            }
            # get the bases that siteB and C are changing into
            baseBC <- paste0(siteB_change[2], siteC_change[2])
            if (baseBC %in% siteBC_coev_pair) {
                Qco3[i, j] <- paste0(Qco3[i, j] , 'f')
            }
            next
        }
    }
}