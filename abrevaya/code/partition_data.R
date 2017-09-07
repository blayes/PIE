library(matrixStats)

panelData <- read.table("../data/panel.txt", header = TRUE)

ncv <- 10

birth <- panelData[panelData$proxy_exists == 1  & (as.numeric(panelData$smoke) == 1), ]
birth$proxy_exists <- NULL
birth$proxy_or_proxyhat <- NULL
birth$stateres <- NULL

momTbl <- table(birth$momid3)
selMom <- names(momTbl)[which(as.numeric(momTbl) == 2)]
selIdx <- as.character(birth$momid3) %in% selMom
birth <- birth[((selIdx & birth$collgrad == 0) & birth$somecoll == 0), ]

birth$novisit <- NULL

saveRDS(birth, "../data/birth.rds")

nameRanefs <- c("dmage", "nlbnl", "gestat", "agesq")
nameFixefs <- c("dmage", "nlbnl", "gestat", "male", "married", "hsgrad", "agesq", "black", "adeqcode2", "adeqcode3", "pretri2", "pretri3")

xmat <- as.matrix(birth[ , nameFixefs])
xmat[ , c("dmage", "nlbnl", "gestat", "agesq")] <- xmat[  , c("dmage", "nlbnl", "gestat", "agesq")] / matrix(colMads(xmat[ , c("dmage", "nlbnl", "gestat", "agesq")]), ncol = 4, nrow = nrow(xmat), byrow = TRUE)

zmat <- as.matrix(birth[ , nameRanefs])
zmat <- zmat / matrix(colMads(zmat), nrow = nrow(zmat), ncol = ncol(zmat), byrow = TRUE)

yvec <- as.numeric(birth[ , "dbirwt"])
yvec <- yvec / mad(yvec)

group <- as.integer(birth[ , "momid3"])

cvidx <- sample(1:ncv, length(unique(group)), replace = TRUE)

unqGrp <- unique(group)
cvdata <- vector("list", ncv)
for (ii in 1:ncv) {
    cvgrps <- unqGrp[cvidx == ii]
    cvdata[[ii]]$x <- cbind("inter" = 1, xmat[group %in% cvgrps, ])
    cvdata[[ii]]$z <- zmat[group %in% cvgrps, ]
    cvdata[[ii]]$group <- group[group %in% cvgrps]
    cvdata[[ii]]$y <- yvec[group %in% cvgrps]
}

train <- vector("list", ncv)
test <- vector("list", ncv)
for (ii in 1:ncv) {
    test[[ii]] <- cvdata[[ii]]
    train[[ii]]$x <- do.call(rbind, lapply(cvdata[-ii], '[[', 'x'))
    train[[ii]]$z <- do.call(rbind, lapply(cvdata[-ii], '[[', 'z'))
    train[[ii]]$group <- do.call(c, lapply(cvdata[-ii], '[[', 'group'))
    train[[ii]]$y <- do.call(c, lapply(cvdata[-ii], '[[', 'y'))
}

saveRDS(train, "../data/train_bwt.rds")
saveRDS(test, "../data/test_bwt.rds")

## wasp partitions
rm(list=ls())

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
