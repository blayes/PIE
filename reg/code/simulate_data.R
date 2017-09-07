rm(list=ls())
set.seed(12345)

setwd("~/pie/reg/code")

repData <- readRDS("/nfsscratch/Users/ssrivastva/pie/reg/data/reg.rds")
nobs <- c(1e4, 1e5)
ndim <- c(10, 100, 200, 300, 400, 500)

for (ii in 1:10) {
  for (jj in 1:6) {
    for (kk in 1:2) {
      partData <- vector("list", 10)
      names(partData) <- paste0("k", 1:10)
      idx <- 1:nrow(repData[[ii]][[jj]][[kk]]$x)
      pidx <- sample(1:10, size = nrow(repData[[ii]][[jj]][[kk]]$x), replace = TRUE)
      for (ll in 1:10) {
        partData[[ll]] <- list(y = repData[[ii]][[jj]][[kk]]$y[pidx == ll],
                               x = repData[[ii]][[jj]][[kk]]$x[pidx == ll, ],
                               mu = repData[[ii]][[jj]][[kk]]$mu[pidx == ll],
                               nrep = nobs[kk] / sum(pidx == ll))
      }
      saveRDS(partData, paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/wasp_reg_cv_", ii, "_p_", ndim[jj], "_n_", nobs[kk], "_k_10", ".rds"))
      cat("ii: ", ii, "jj: ", jj, "kk: ", kk, "\n")
    }
  }
}

rm(list=ls())
set.seed(12345)

setwd("~/pie/reg/code")

repData <- readRDS("/nfsscratch/Users/ssrivastva/pie/reg/data/reg.rds")
nobs <- c(1e4, 1e5)
ndim <- c(10, 100, 200, 300, 400, 500)

for (ii in 1:10) {
  for (jj in 1:6) {
    for (kk in 1:2) {
      partData <- vector("list", 20)
      names(partData) <- paste0("k", 1:20)
      idx <- 1:nrow(repData[[ii]][[jj]][[kk]]$x)
      pidx <- sample(1:20, size = nrow(repData[[ii]][[jj]][[kk]]$x), replace = TRUE)
      for (ll in 1:20) {
        partData[[ll]] <- list(y = repData[[ii]][[jj]][[kk]]$y[pidx == ll],
                               x = repData[[ii]][[jj]][[kk]]$x[pidx == ll, ],
                               mu = repData[[ii]][[jj]][[kk]]$mu[pidx == ll],
                               nrep = nobs[kk] / sum(pidx == ll))
      }
      saveRDS(partData, paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/wasp_reg_cv_", ii, "_p_", ndim[jj], "_n_", nobs[kk], "_k_20", ".rds"))
      cat("ii: ", ii, "jj: ", jj, "kk: ", kk, "\n")
    }
  }
}
