cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

setwd('/Users/ssrivastva/pie/abrevaya/code/')

if (mtd == 1) {
    source("ml_estimation.R")
    cvtrain <- readRDS("../data/train_bwt.rds")
    train <- cvtrain[[id]]
    res <- estimateMLE(train$y, train$x, train$z, train$group, family = "gaussian")
    fname <- paste0("/nfsscratch/Users/ssrivastva/pie/abrevaya/result/ml/ml_res_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("mcmc_sampler.R")
    cvtrain <- readRDS("../data/train_bwt.rds")
    train <- cvtrain[[id]]
    res <- sampleFromMixMdl(train$y, train$x, train$z, train$group, 10000, 5000, 5, dirp = '/nfsscratch/Users/ssrivastva/pie/abrevaya/', id)
  } else if (mtd == 3) {
    source("variational_bayes.R")
    cvtrain <- readRDS("../data/train_bwt.rds")
    train <- cvtrain[[id]]
    res <- fitLinearMixefEffectsVB(train$y, train$x, train$z, train$group, 100)
    fname <- paste0("/nfsscratch/Users/ssrivastva/pie/abrevaya/result/vb/vb_res_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 4) {
    source("wasp_sampler.R")
    reps <- rep(1:10, each = 20)
    subs <- rep(1:20, times = 10)
    cvtrain <- readRDS("../data/wasp_train_bwt.rds")
    train <- cvtrain[[reps[id]]][[subs[id]]]
    res <- sampleFromWaspMixMdl(train$y, train$x, train$z, train$group, 20, 10000, 5000, 5, dirp = '/nfsscratch/Users/ssrivastva/pie/abrevaya/', reps[id], subs[id])
} else if (mtd == 5) {
    source("comp_sampler.R")

    reps <- rep(1:10, each = 20)
    subs <- rep(1:20, times = 10)

    cvtrain <- readRDS("../data/wasp_train_bwt.rds")
    train <- cvtrain[[reps[id]]][[subs[id]]]

    res <- sampleFromCompMixMdl(train$y, train$x, train$z, train$group, 20, 10000, 5000, 5, id)
    fname <- paste0("/Shared/ssrivastva/pie/abrevaya/result/comp/samp/cv_", reps[id], "_k_", subs[id], ".rds")
    saveRDS(res, fname)
} else (mtd == 6) {
    library(parallelMCMCcombine)

    cid <- id

    subfix <- array(0.0, dim = c(14, 1000, 20))
    subran <- array(0.0, dim = c(6, 1000, 20))
    tmp <- numeric(20)

    for (kk in 1:20) {
        fname <- paste0("/Shared/ssrivastva/pie/abrevaya/result/comp/samp/cv_", cid, "_k_", kk, ".rds")
        samp <- readRDS(fname)
        cnames <- colnames(samp$samples)
        subfix[ , , kk] <- t(samp$samples[ , 1:14])
        tmpvar <- samp$samples[ ,  grep("ranef", cnames, ignore.case=T)[c(1:3, 5, 6, 9)]]
        subran[ , , kk] <- t(tmpvar / cbind(1, sqrt(tmpvar[ , 1] * tmpvar[ , 4]), sqrt(tmpvar[ , 1] * tmpvar[ , 6]), 1, sqrt(tmpvar[ , 4] * tmpvar[ , 6]), 1))
        tmp[kk] <- samp$time[3]
    }

    strt1 <- proc.time()
    scottFix <- consensusMCindep(subchain = subfix)
    scottRan <- consensusMCindep(subchain = subran)
    end1 <- proc.time()
    strt2 <- proc.time()
    xingFix <- semiparamDPE(subchain = subfix)
    xingRan <- semiparamDPE(subchain = subran)
    end2 <- proc.time()

    fname1 <- paste0("/Shared/ssrivastva/pie/abrevaya/result/cons/corr_cv_", cid, "_k20.rds")
    fname2 <- paste0("/Shared/ssrivastva/pie/abrevaya/result/xing/corr_cv_", cid, "_k20.rds")

    saveRDS(list(fix = t(scottFix), ran = t(scottRan), time = mean(tmp) + end1[3] - strt1[3]), fname1)
    saveRDS(list(fix = t(xingFix), ran = t(xingRan), time = mean(tmp) + end2[3] - strt2[3]), fname2)
}
