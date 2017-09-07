cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
  cvs <- rep(1:10, each = 12)
  ndims <- rep(rep(1:6, each = 2), times = 10)
  nobss <- rep(rep(1:2, times = 6), times = 10)

  cid <- cvs[id]
  did <- ndims[id]
  nid <- nobss[id]

  cvtrain <- readRDS("/nfsscratch/Users/ssrivastva/pie/reg/data/reg.rds")
  train <- cvtrain[[cid]][[did]][[nid]]
  rm(cvtrain)

  startTime <- proc.time()
  res <- lm(y ~ 0 + x, data = train)
  endTime <- proc.time()

  fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/result/lap/lap_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], ".rds")
  saveRDS(list(res = res, time = endTime - startTime), fname)
} else if (mtd == 2) {
  source("mcmc_sampler.R")

  cvs <- rep(1:10, each = 12)
  ndims <- rep(rep(1:6, each = 2), times = 10)
  nobss <- rep(rep(1:2, times = 6), times = 10)

  cid <- cvs[id]
  did <- ndims[id]
  nid <- nobss[id]

  cvtrain <- readRDS("/nfsscratch/Users/ssrivastva/pie/reg/data/reg.rds")
  train <- cvtrain[[cid]][[did]][[nid]]
  rm(cvtrain)

  res <- sampleGdp(train$y, train$x, 100, 10000, 5000, 5)
  fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/result/mcmc/mcmc_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], ".rds")
  saveRDS(res, fname)
} else if (mtd == 3) {
  ## sub = 10
  cvs <- rep(1:10, each = 12)
  ndims <- rep(rep(1:6, each = 2), times = 10)
  nobss <- rep(rep(1:2, times = 6), times = 10)

  tmp <- cbind(cvs, ndims, nobss)
  wids <- cbind(tmp[rep(1:nrow(tmp), each = 10), ], subs = rep(1:10, times = 120))

  cid <- wids[id, 1]
  did <- wids[id, 2]
  nid <- wids[id, 3]
  sid <- wids[id, 4]

  cvtrain <- readRDS(paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/wasp_reg_cv_", cid,"_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_10", ".rds"))

  train <- cvtrain[[sid]]
  rm(cvtrain)

  startTime <- proc.time()
  res <- lm(y ~ 0 + x, data = train)
  endTime <- proc.time()
  fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/samp/k10/", cid, "/lap_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_", sid, ".rds")
  saveRDS(list(res = res, time = endTime - startTime), fname)
} else if (mtd == 4) {
  source("wasp_sampler.R")
  cvs <- rep(1:10, each = 12)
  ndims <- rep(rep(1:6, each = 2), times = 10)
  nobss <- rep(rep(1:2, times = 6), times = 10)

  tmp <- cbind(cvs, ndims, nobss)
  wids <- cbind(tmp[rep(1:nrow(tmp), each = 10), ], rep(1:10, times = 120))

  cid <- wids[id, 1]
  did <- wids[id, 2]
  nid <- wids[id, 3]
  sid <- wids[id, 4]

  cvtrain <- readRDS(paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/wasp_reg_cv_", cid,"_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_10", ".rds"))

  train <- cvtrain[[sid]]
  rm(cvtrain)

  res <- sampleWaspGdp(train$y, train$x, train$nrep, 100, 10000, 5000, 5)
  fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/samp/k10/", cid, "/wasp_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_", sid, ".rds")
  saveRDS(res, fname)
} else if (mtd == 5) {
  source("wasp_sampler.R")
  cvs <- rep(1:10, each = 12)
  ndims <- rep(rep(1:6, each = 2), times = 10)
  nobss <- rep(rep(1:2, times = 6), times = 10)

  tmp <- cbind(cvs, ndims, nobss)
  wids <- cbind(tmp[rep(1:nrow(tmp), each = 20), ], rep(1:20, times = 120))
  cid <- wids[id, 1]
  did <- wids[id, 2]
  nid <- wids[id, 3]
  sid <- wids[id, 4]

  cvtrain <- readRDS(paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/wasp_reg_cv_", cid,"_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_20", ".rds"))

  train <- cvtrain[[sid]]
  rm(cvtrain)

  res <- sampleWaspGdp(train$y, train$x, train$nrep, 100, 10000, 5000, 5)
  fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/samp/k20/", cid, "/wasp_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_", sid, ".rds")
  saveRDS(res, fname)
} else if (mtd == 6) {
  ## sub = 20
  cvs <- rep(1:10, each = 12)
  ndims <- rep(rep(1:6, each = 2), times = 10)
  nobss <- rep(rep(1:2, times = 6), times = 10)

  tmp <- cbind(cvs, ndims, nobss)
  wids <- cbind(tmp[rep(1:nrow(tmp), each = 20), ], subs = rep(1:20, times = 120))
  cid <- wids[id, 1]
  did <- wids[id, 2]
  nid <- wids[id, 3]
  sid <- wids[id, 4]

  cvtrain <- readRDS(paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/wasp_reg_cv_", cid,"_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_20", ".rds"))

  train <- cvtrain[[sid]]
  rm(cvtrain)

  startTime <- proc.time()
  res <- lm(y ~ 0 + x, data = train)
  endTime <- proc.time()
  fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/samp/k20/", cid, "/lap_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_", sid, ".rds")
  saveRDS(list(res = res, time = endTime - startTime), fname)
} else if (mtd == 7) {
  source("obtain_barycenter.R")
  ## sub = 10
  cvs <- rep(1:10, each = 12)
  ndims <- rep(rep(1:6, each = 2), times = 10)
  nobss <- rep(rep(1:2, times = 6), times = 10)

  wids <- cbind(cvs, ndims, nobss)

  cid <- wids[id, 1]
  did <- wids[id, 2]
  nid <- wids[id, 3]

  muList <- vector("list", 10)
  sigList <- vector("list", 10)
  rtime <- rep(0, 10)
  for (sid in 1:10) {
    fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/samp/k10/", cid, "/lap_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_", sid, ".rds")
    dat <- readRDS(fname)
    muList[[sid]] <- coef(dat$res)
    sigList[[sid]] <- vcov(dat$res) / 10
    rtime[sid] <- dat$time[3]
  }
  startTime <- proc.time()
  res <- obtainGaussianBarycenter(muList, sigList, wts = rep(1/10, 10))
  endTime <- proc.time()

  fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/result/ml/k10/lap_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k10.rds")
  saveRDS(list(res = res, time = mean(rtime) + endTime[3] - startTime[3]), fname)
} else {
  source("obtain_barycenter.R")
  ## sub = 20
  cvs <- rep(1:10, each = 12)
  ndims <- rep(rep(1:6, each = 2), times = 10)
  nobss <- rep(rep(1:2, times = 6), times = 10)

  wids <- cbind(cvs, ndims, nobss)

  cid <- wids[id, 1]
  did <- wids[id, 2]
  nid <- wids[id, 3]

  muList <- vector("list", 20)
  sigList <- vector("list", 20)
  rtime <- rep(0, 20)
  for (sid in 1:20) {
    fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/data/samp/k20/", cid, "/lap_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k_", sid, ".rds")
    dat <- readRDS(fname)
    muList[[sid]] <- coef(dat$res)
    sigList[[sid]] <- vcov(dat$res) / 20
    rtime[sid] <- dat$time[3]
  }
  startTime <- proc.time()
  res <- obtainGaussianBarycenter(muList, sigList, wts = rep(1/20, 20))
  endTime <- proc.time()

  fname <- paste0("/nfsscratch/Users/ssrivastva/pie/reg/result/ml/k20/lap_cv_", cid, "_p_", c(10, 100, 200, 300, 400, 500)[did], "_n_", c(1e4, 1e5)[nid], "_k20.rds")
  saveRDS(list(res = res, time = mean(rtime) + endTime[3] - startTime[3]), fname)
}
