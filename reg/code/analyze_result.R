rm(list=ls())

setwd("/nfsscratch/Users/ssrivastva/pie/reg/result/")

library(KernSmooth)
library(MCMCpack)
library(matrixStats)
library(RColorBrewer)
colors <- brewer.pal(6, "Set1")

mcmcCoef <- list()
lapCoef <- list()
mlCoef <- list()
pieCoef <- list()

mcmcTime <- list()
lapTime <- list()
mlTime <- list()
pieTime <- list()

nobs <- c(10000, 1e05)
ndim <- c(10, 100, 200, 300, 400, 500)
nsubs <- c(10, 20)

for (cc in 1:10) {
    mcmcCoef[[cc]] <- list()
    mcmcTime[[cc]] <- list()
    for (pp in 1:6) {
        mcmcCoef[[cc]][[pp]] <- list()
        mcmcTime[[cc]][[pp]] <- list()
        for (nn in 1:2) {
            dat <- readRDS(paste0("mcmc/mcmc_cv_", cc, "_p_", ndim[pp], "_n_", nobs[nn], ".rds"))
            mcmcCoef[[cc]][[pp]][[nn]] <- dat$betaSamp
            mcmcTime[[cc]][[pp]][[nn]] <- dat$time[3]
        }
    }
}

for (cc in 1:10) {
    mlCoef[[cc]] <- list()
    mlTime[[cc]] <- list()
    for (pp in 1:6) {
        mlCoef[[cc]][[pp]] <- list()
        mlTime[[cc]][[pp]] <- list()
        for (nn in 1:2) {
            if (nn == 1 & pp == 6) {
                dat1 <- readRDS(paste0("ml/k10/lap_cv_", cc, "_p_", ndim[pp], "_n_", nobs[nn], "_k10.rds"))
                bb1 <- dat1$res$mean
                vv1 <- dat1$res$cov
                mlCoef[[cc]][[pp]][[nn]] <- list(t(crossprod(chol(vv1), matrix(rnorm(1000 * nrow(vv1)), nrow(vv1), 1000)) + bb1),
                                             t(crossprod(chol(vv1), matrix(rnorm(1000 * nrow(vv1)), nrow(vv1), 1000)) + bb1)
                                             )
                mlTime[[cc]][[pp]][[nn]] <- list(dat1$time[1], dat1$time[1])
            } else {
                dat1 <- readRDS(paste0("ml/k10/lap_cv_", cc, "_p_", ndim[pp], "_n_", nobs[nn], "_k10.rds"))
                bb1 <- dat1$res$mean
                vv1 <- dat1$res$cov
                dat2 <- readRDS(paste0("ml/k20/lap_cv_", cc, "_p_", ndim[pp], "_n_", nobs[nn], "_k20.rds"))
                bb2 <- dat2$res$mean
                vv2 <- dat2$res$cov
                mlCoef[[cc]][[pp]][[nn]] <- list(t(crossprod(chol(vv1), matrix(rnorm(1000 * nrow(vv1)), nrow(vv1), 1000)) + bb1),
                                             t(crossprod(chol(vv2), matrix(rnorm(1000 * nrow(vv2)), nrow(vv2), 1000)) + bb2)
                                             )
                mlTime[[cc]][[pp]][[nn]] <- list(dat1$time[1], dat2$time[1])
            }
        }
    }
}

for (cc in 1:10) {
    lapCoef[[cc]] <- list()
    lapTime[[cc]] <- list()
    for (pp in 1:6) {
        lapCoef[[cc]][[pp]] <- list()
        lapTime[[cc]][[pp]] <- list()
        for (nn in 1:2) {
            dat <- readRDS(paste0("lap/lap_cv_", cc, "_p_", ndim[pp], "_n_", nobs[nn], ".rds"))
            bb <- coef(dat$res)
            vv <- vcov(dat$res)
            lapCoef[[cc]][[pp]][[nn]] <- t(crossprod(chol(vv), matrix(rnorm(1000 * nrow(vv)), nrow(vv), 1000)) + bb)
            lapTime[[cc]][[pp]][[nn]] <- dat$time[3]
        }
    }
}

for (cc in 1:10) {
    pieCoef[[cc]] <- list()
    pieTime[[cc]] <- list()
    for (pp in 1:6) {
        pieCoef[[cc]][[pp]] <- list()
        pieTime[[cc]][[pp]] <- list()
        for (nn in 1:2) {
            pieCoef[[cc]][[pp]][[nn]] <- list()
            pieTime[[cc]][[pp]][[nn]] <- list()
            for (kk in 1:2) {
                pieCoef[[cc]][[pp]][[nn]][[kk]] <- list()
                tmp <- numeric(0)
                for (ss in 1:nsubs[kk]) {
                    dat <- readRDS(paste0("../data/samp/k", nsubs[kk], "/", cc, "/wasp_cv_", cc, "_p_", ndim[pp], "_n_", nobs[nn], "_k_", ss, ".rds"))
                    pieCoef[[cc]][[pp]][[nn]][[kk]][[ss]] <- dat$betaSamp
                    tmp <- append(tmp, dat$time[3])
                }
                pieTime[[cc]][[pp]][[nn]][[kk]] <- mean(tmp)
            }
        }
    }
}

pieCoefEst <- list()
for (cc in 1:10) {
    pieCoefEst[[cc]] <- list()
    for (pp in 1:6) {
        pieCoefEst[[cc]][[pp]] <- list()
        for (nn in 1:2) {
            pieCoefEst[[cc]][[pp]][[nn]] <- list()
            for (kk in 1:2) {
                pieCoefEst[[cc]][[pp]][[nn]][[kk]] <- list()
                tmp <- list()
                for (ss in 1:nsubs[kk]) {
                    tmp[[ss]] <- do.call(cbind, lapply(split(pieCoef[[cc]][[pp]][[nn]][[kk]][[ss]], col(pieCoef[[cc]][[pp]][[nn]][[kk]][[ss]])),
                                                       function(x) quantile(x, prob = seq(0.0, 1, length = 1000))))
                }
                pieCoefEst[[cc]][[pp]][[nn]][[kk]] <- matrix(0.0, nrow = dim(tmp[[1]])[1], ncol = dim(tmp[[1]])[2])
                for (dd in 1:ncol(pieCoefEst[[cc]][[pp]][[nn]][[kk]])) {
                    pieCoefEst[[cc]][[pp]][[nn]][[kk]][ , dd] <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[ , dd])))
                }
            }
        }
    }
}

################################## plots and tables ############################

resCoef <- vector("list", 4)
names(resCoef) <- c("MCMC", "ML", "PIE", "Laplace")
for (cc in 1:10) {
    resCoef[["MCMC"]][[cc]] <- list()
    resCoef[["ML"]][[cc]] <- list()
    resCoef[["PIE"]][[cc]] <- list()
    resCoef[["Laplace"]][[cc]] <- list()
    for (pp in 1:6) {
        resCoef[["MCMC"]][[cc]][[pp]] <- list()
        resCoef[["ML"]][[cc]][[pp]] <- list()
        resCoef[["PIE"]][[cc]][[pp]] <- list()
        resCoef[["Laplace"]][[cc]][[pp]] <- list()
        for (nn in 1:2) {
            resCoef[["MCMC"]][[cc]][[pp]][[nn]] <- list()
            resCoef[["ML"]][[cc]][[pp]][[nn]] <- list()
            resCoef[["PIE"]][[cc]][[pp]][[nn]] <- list()
            resCoef[["Laplace"]][[cc]][[pp]][[nn]] <- list()
            for (kk in 1:2) {
                resCoef[["MCMC"]][[cc]][[pp]][[nn]][[kk]] <- list()
                resCoef[["ML"]][[cc]][[pp]][[nn]][[kk]] <- list()
                resCoef[["PIE"]][[cc]][[pp]][[nn]][[kk]] <- list()
                resCoef[["Laplace"]][[cc]][[pp]][[nn]][[kk]] <- list()
                for (jj in 1:ncol(mcmcCoef[[cc]][[pp]][[nn]])) {
                    rr <- range(c(mcmcCoef[[cc]][[pp]][[nn]][ , jj],
                                  mlCoef[[cc]][[pp]][[nn]][[kk]][ , jj],
                                  pieCoefEst[[cc]][[pp]][[nn]][[kk]][ , jj],
                                  lapCoef[[cc]][[pp]][[nn]][ , jj]))
                    bw1 <- dpik(mcmcCoef[[cc]][[pp]][[nn]][ , jj], range.x = rr)
                    bw2 <- dpik(mlCoef[[cc]][[pp]][[nn]][[kk]][ , jj], range.x = rr)
                    bw3 <- dpik(pieCoefEst[[cc]][[pp]][[nn]][[kk]][ , jj], range.x = rr)
                    bw5 <- dpik(lapCoef[[cc]][[pp]][[nn]][ , jj], range.x = rr)
                    dens1 <- bkde(mcmcCoef[[cc]][[pp]][[nn]][ , jj], bandwidth = bw1, range.x = rr)
                    dens2 <- bkde(mlCoef[[cc]][[pp]][[nn]][[kk]][ , jj], bandwidth = bw2, range.x = rr)
                    dens3 <- bkde(pieCoefEst[[cc]][[pp]][[nn]][[kk]][ , jj], bandwidth = bw3, range.x = rr)
                    dens5 <- bkde(lapCoef[[cc]][[pp]][[nn]][ , jj], bandwidth = bw5, range.x = rr)
                    resCoef[["MCMC"]][[cc]][[pp]][[nn]][[kk]][[jj]] <- dens1
                    resCoef[["ML"]][[cc]][[pp]][[nn]][[kk]][[jj]] <- dens2
                    resCoef[["PIE"]][[cc]][[pp]][[nn]][[kk]][[jj]] <- dens3
                    resCoef[["Laplace"]][[cc]][[pp]][[nn]][[kk]][[jj]] <- dens5
                }
            }
        }
    }
}

accCoef <- list()
for (pp in 1:6) {
    accCoef[[pp]] <- list()
    for (nn in 1:2) {
        accCoef[[pp]][[nn]] <- list()
        for (kk in 1:2) {
            accCoef[[pp]][[nn]][[kk]] <- matrix(0.0, 3, ncol(mcmcCoef[[1]][[pp]][[nn]]))
            for (jj in 1:ncol(mcmcCoef[[1]][[pp]][[nn]])) {
                for (cc in 1:10) {
                    accCoef[[pp]][[nn]][[kk]][1, jj] <- accCoef[[pp]][[nn]][[kk]][1, jj] +
                        (1 - sum(abs(resCoef[["MCMC"]][[cc]][[pp]][[nn]][[kk]][[jj]]$y  - resCoef[["Laplace"]][[cc]][[pp]][[nn]][[kk]][[jj]]$y) * diff(resCoef[["MCMC"]][[cc]][[pp]][[nn]][[kk]][[jj]]$x)[1]) / 2)
                    accCoef[[pp]][[nn]][[kk]][2, jj] <- accCoef[[pp]][[nn]][[kk]][2, jj] +
                        (1 - sum(abs(resCoef[["MCMC"]][[cc]][[pp]][[nn]][[kk]][[jj]]$y  - resCoef[["PIE"]][[cc]][[pp]][[nn]][[kk]][[jj]]$y) * diff(resCoef[["MCMC"]][[cc]][[pp]][[nn]][[kk]][[jj]]$x)[1]) / 2)
                    accCoef[[pp]][[nn]][[kk]][3, jj] <- accCoef[[pp]][[nn]][[kk]][3, jj] +
                        (1 - sum(abs(resCoef[["MCMC"]][[cc]][[pp]][[nn]][[kk]][[jj]]$y  - resCoef[["ML"]][[cc]][[pp]][[nn]][[kk]][[jj]]$y) * diff(resCoef[["MCMC"]][[cc]][[pp]][[nn]][[kk]][[jj]]$x)[1]) / 2)
                }
            }
            accCoef[[pp]][[nn]][[kk]] <- accCoef[[pp]][[nn]][[kk]] / 10
        }
    }
}

accMdl <- list()
for (pp in 1:6) {
    accMdl[[pp]] <- list()
    for (nn in 1:2) {
        accMdl[[pp]][[nn]] <- list()
        zeros1 <- rowMeans(accCoef[[pp]][[nn]][[1]][ , -(1:(0.1 * ncol(accCoef[[pp]][[nn]][[1]]))), drop = FALSE])
        nzeros1 <- rowMeans(accCoef[[pp]][[nn]][[1]][ , (1:(0.1 * ncol(accCoef[[pp]][[nn]][[1]]))), drop = FALSE])
        zeros2 <- rowMeans(accCoef[[pp]][[nn]][[2]][ , -(1:(0.1 * ncol(accCoef[[pp]][[nn]][[2]]))), drop = FALSE])
        nzeros2 <- rowMeans(accCoef[[pp]][[nn]][[2]][ , (1:(0.1 * ncol(accCoef[[pp]][[nn]][[2]]))), drop = FALSE])
        accMdl[[pp]][[nn]] <- round(rbind(cbind(zero = zeros1, nzero = nzeros1),
                                          cbind(zero = zeros2, nzero = nzeros2)[-1, ]), 2)
        rownames(accMdl[[pp]][[nn]]) <- c("Laplace", "PIE(k=10)", "ML-II(k=10)", "PIE(k=20)", "ML-II(k=20)")
    }
    names(accMdl[[pp]]) <- c("n=10k", "n=100k")
}
names(accMdl) <- c("p=10", "p=100", "p=200", "p=300", "p=400", "p=500")
