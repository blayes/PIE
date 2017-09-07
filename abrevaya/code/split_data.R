## wasp partitions
rm(list=ls())

set.seed(54321)

npart <- 20
nrep <- 10

reps <- readRDS("../data/train_bwt.rds")
parts <- vector("list", length = nrep)
partsIdx <- vector("list", length = nrep)
for (ii in 1:nrep) {
    parts[[ii]] <- vector("list", length = npart)
}

for (r in seq_len(nrep)) {
    unqGrp <- unique(reps[[r]]$group)
    pidx <- sample(1:npart, length(unqGrp), replace = TRUE)
    grp <- reps[[r]]$group
    for (ii in 1:npart) {
        pgrps <- unqGrp[pidx == ii]
        idx <- grp %in% pgrps
        parts[[r]][[ii]]$nobs <- sum(idx)
        parts[[r]][[ii]]$group <- reps[[r]]$group[idx]
        parts[[r]][[ii]]$x <- reps[[r]]$x[idx, ]
        parts[[r]][[ii]]$y <- reps[[r]]$y[idx]
        parts[[r]][[ii]]$z <- reps[[r]]$z[idx, ]
    }
}

saveRDS(parts, "../data/wasp_train_bwt.rds")
