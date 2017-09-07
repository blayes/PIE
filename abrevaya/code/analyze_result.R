rm(list=ls())

setwd("/Shared/ssrivastva/pie/abrevaya/result/")

library(KernSmooth)
library(MCMCpack)
library(matrixStats)
library(RColorBrewer)
colors <- brewer.pal(6, "Set1")

mcmcFix <- list()
mlFix <- list()
vbFix <- list()
waspFix <- list()
pieFix <- list()
xingFix <- list()
scotFix <- list()

mcmcRan <- list()
mlRan <- list()
vbRan <- list()
waspRan <- list()
pieRan <- list()
xingRan <- list()
scotRan <- list()

mcmcCorr <- list()
mlCorr <- list()
vbCorr <- list()
waspCorr <- list()
pieCorr <- list()
xingCorr <- list()
scotCorr <- list()

mcmcTime <- list()
mlTime <- list()
vbTime <- list()
waspTime <- list()
pieTime <- list()
xingTime <- list()
scotTime <- list()

for (ii in 1:10) {
    dat <- readRDS(paste0("xing/cv_", ii, "_k20.rds"))
    xingFix[[ii]] <- dat$fix
    xingRan[[ii]] <- dat$ran
    xingTime[[ii]] <- dat$time
}

for (ii in 1:10) {
    dat <- readRDS(paste0("cons/cv_", ii, "_k20.rds"))
    scotFix[[ii]] <- dat$fix
    scotRan[[ii]] <- dat$ran
    scotTime[[ii]] <- dat$time
}

for (ii in 1:10) {
    dat <- readRDS(paste0("xing/corr_cv_", ii, "_k20.rds"))
    xingCorr[[ii]] <- dat$ran
}

for (ii in 1:10) {
    dat <- readRDS(paste0("cons/corr_cv_", ii, "_k20.rds"))
    scotCorr[[ii]] <- dat$ran
}

for (ii in 1:10) {
    dat <- readRDS(paste0("mcmc/samp_time_", ii, ".rds"))
    mcmcFix[[ii]] <- dat$samples[ , 1:14]
    mcmcRan[[ii]] <- dat$samples[ , c(15:17, 19:20, 23)]
    mcmcTime[[ii]] <- dat$time[3]
}

for (ii in 1:10) {
    dat <- readRDS(paste0("vb/vb_res_", ii, ".rds"))
    vbFix[[ii]] <- t(crossprod(chol(dat$coef$cov), matrix(rnorm(nrow(dat$coef$cov) * 1000), nrow(dat$coef$cov), 1000)) + dat$coef$mu)
    tmp <- array(0.0, dim = c(dim(dat$cov$scale), 1000))
    for (ll in 1:1000) {
        tmp[ , , ll] <- riwish(dat$cov$df, dat$cov$scale)
    }
    vbRan[[ii]] <- t(apply(tmp, 3,
                           function(x) c(x[1, 1], x[2, 1], x[3, 1],
                                         x[2, 2], x[2, 3],
                                         x[3, 3])))
    vbTime[[ii]] <- dat$time[3]
}

for (ii in 1:10) {
    dat <- readRDS(paste0("../result/ml/ml_res_", ii, ".rds"))
    mlFix[[ii]] <- t(crossprod(chol(dat$fixefCov), matrix(rnorm(nrow(dat$fixefCov) * 1000), nrow(dat$fixefCov), 1000)) + dat$fixef)
    mlRan[[ii]] <- c(dat$ranefCov[1, 1], dat$ranefCov[1, 2], dat$ranefCov[1, 3],
                     dat$ranefCov[2, 2], dat$ranefCov[2, 3], dat$ranefCov[3, 3])
    mlCorr[[ii]] <- c(dat$ranefCov[1, 1],
                      dat$ranefCov[1, 2] / sqrt(dat$ranefCov[1, 1] * dat$ranefCov[2, 2]),
                      dat$ranefCov[1, 3] / sqrt(dat$ranefCov[1, 1] * dat$ranefCov[3, 3]),
                      dat$ranefCov[2, 2], dat$ranefCov[2, 3],
                      dat$ranefCov[3, 3] / sqrt(dat$ranefCov[2, 2] * dat$ranefCov[3, 3]))
    mlTime[[ii]] <- dat$time[3]
}

for (ii in 1:10) {
    waspFix[[ii]] <- matrix(0.0, 1000, 14)
    for (jj in 1:14) {
        dat <- read.csv(paste0("wasp/fixef_", jj, "_rep_", ii, ".csv"), header = FALSE)
        pp <- as.numeric(dat[ , 2])
        pp[pp < 1e-10] <- 0.0
        waspFix[[ii]][ , jj] <- sample(as.numeric(dat[ , 1]), 1000, replace = TRUE, prob = pp)
    }

    waspRan[[ii]] <- matrix(0.0, 1000, 6)
    for (jj in 1:6) {
        dat <- read.csv(paste0("wasp/ranef_", jj, "_rep_", ii, ".csv"), header = FALSE)
        pp <- as.numeric(dat[ , 2])
        pp[pp < 1e-10] <- 0.0
        waspRan[[ii]][ , jj] <- sample(as.numeric(dat[ , 1]), 1000, replace = TRUE, prob = pp)
    }

    waspCorr[[ii]] <- matrix(NA, 1000, 6)
    for (jj in c(1, 4, 6)) {
        dat <- read.csv(paste0("wasp/ranef_", jj, "_rep_", ii, ".csv"), header = FALSE)
        pp <- as.numeric(dat[ , 2])
        pp[pp < 1e-10] <- 0.0
        waspCorr[[ii]][ , jj] <- sample(as.numeric(dat[ , 1]), 1000, replace = TRUE, prob = pp)
    }

    dat1 <- read.csv(paste0("wasp/ranef_times_rep_", ii, ".csv"), header = FALSE)
    dat2 <- read.csv(paste0("wasp/fixef_times_rep_", ii, ".csv"), header = FALSE)

    tmp <- rep(0, 20)
    for (jj in 1:20) {
        dat3 <- readRDS(paste0("../data/wasp_samp/", ii, "/samp_time_", jj, ".rds"))
        tmp[jj] <- dat3$time[3]
    }
    waspTime[[ii]] <- mean(c(unlist(dat1), unlist(dat2))) + mean(tmp)
}

for (ii in 1:10) {
    pieFix[[ii]] <- list()
    pieRan[[ii]] <- list()
    pieCorr[[ii]] <- list()
    tmp <- rep(0, 20)
    for (jj in 1:20) {
        dat <- readRDS(paste0("../data/wasp_samp/", ii, "/samp_time_", jj, ".rds"))
        pieFix[[ii]][[jj]] <- dat$samples[ , 1:14]
        pieRan[[ii]][[jj]] <- dat$samples[ , c(15:17, 19:20, 23)]
        pieCorr[[ii]][[jj]] <- pieRan[[ii]][[jj]] / cbind(1, sqrt(pieRan[[ii]][[jj]][ , 1] * pieRan[[ii]][[jj]][ , 4]), sqrt(pieRan[[ii]][[jj]][ , 1] * pieRan[[ii]][[jj]][ , 6]), 1, sqrt(pieRan[[ii]][[jj]][ , 4] * pieRan[[ii]][[jj]][ , 6]), 1)
        tmp[jj] <- dat$time[3]
    }
    pieTime[[ii]] <- mean(tmp)
}

pieRanEst <- list()
pieCorrEst <- list()
pieFixEst <- list()
for (ii in 1:10) {
    tmp <- list()
    for (jj in 1:20) {
        tmp[[jj]] <- do.call(cbind, lapply(split(pieRan[[ii]][[jj]], col(pieRan[[ii]][[jj]])),
                                           function(x) quantile(x, prob = seq(0.0, 1, length = 1000))))
    }

    pieRanEst[[ii]] <- matrix(0.0, nrow = dim(tmp[[1]])[1], ncol = dim(tmp[[1]])[2])
    for (jj in 1:ncol(pieRanEst[[ii]])) {
        pieRanEst[[ii]][ , jj] <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[ , jj])))
    }

    tmp <- list()
    for (jj in 1:20) {
        tmp[[jj]] <- do.call(cbind, lapply(split(pieCorr[[ii]][[jj]], col(pieCorr[[ii]][[jj]])),
                                           function(x) quantile(x, prob = seq(0.0, 1, length = 1000))))
    }

    pieCorrEst[[ii]] <- matrix(0.0, nrow = dim(tmp[[1]])[1], ncol = dim(tmp[[1]])[2])
    for (jj in 1:ncol(pieCorrEst[[ii]])) {
        pieCorrEst[[ii]][ , jj] <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[ , jj])))
    }

    tmp <- list()
    for (jj in 1:20) {
        tmp[[jj]] <- do.call(cbind, lapply(split(pieFix[[ii]][[jj]], col(pieFix[[ii]][[jj]])),
                                           function(x) quantile(x, prob = seq(0.0, 1, length = 1000))))
    }

    pieFixEst[[ii]] <- matrix(0.0, nrow = dim(tmp[[1]])[1], ncol = dim(tmp[[1]])[2])
    for (jj in 1:ncol(pieFixEst[[ii]])) {
        pieFixEst[[ii]][ , jj] <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[ , jj])))
    }
}

mcmcCorr <- lapply(mcmcRan, function(x) {
    x / cbind(1, sqrt(x[ , 1] * x[ , 4]), sqrt(x[ , 1] * x[ , 6]), 1, sqrt(x[ , 4] * x[ , 6]), 1)
})

vbCorr <- lapply(vbRan, function(x) {
    x / cbind(1, sqrt(x[ , 1] * x[ , 4]), sqrt(x[ , 1] * x[ , 6]), 1, sqrt(x[ , 4] * x[ , 6]), 1)
})

dims <- c(NA, 1, 2, NA, 3, NA)
for (cc in 1:10) {
    for (dd in c(2, 3, 5)) {
        cat("loaded: ", paste0("wasp/corr_", dims[dd], "_rep_", cc, ".csv"), "\n")
        dat <- read.csv(paste0("wasp/corr_", dims[dd], "_rep_", cc, ".csv"), header = FALSE)
        pp <- as.numeric(dat[ , 2])
        pp[pp < 1e-10] <- 0.0
        waspCorr[[cc]][ , dd] <- sample(as.numeric(dat[ , 1]), 1000, replace = TRUE, prob = pp)
    }
}

############################# plots and tables ##################################
rnames <- c("(dmage, dmage)",  "(nlbnl, nlbnl)", "(gestat, gestat)",
            "(dmage, nbnl)",  "(dmage, gestat)", "(nlbnl, gestat)")
idxran <- c(1, 4, 6, 2, 3, 5)

idxrf <- 1
datRan <- list(MCMC = mcmcRan[[idxrf]][ , idxran],
               VB = vbRan[[idxrf]][ , idxran],
               CMC = scotRan[[idxrf]][ , idxran],
               SC = xingRan[[idxrf]][ , idxran],
               WASP = waspRan[[idxrf]][ , idxran],
               PIE = pieRanEst[[idxrf]][ , idxran]
               )

pdf("~/pie/abrevaya/result/img/abe_ran_box.pdf", 30, 10)
par(mfrow = c(2, 3))
par(cex = 1)
par(mar = c(2, 5, 3, 0.5), oma = c(1, 0, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (dd in 1:length(idxran)) {
    rr <- range(unlist(lapply(datRan, function(x) x[ , dd]))) + c(-0.001, 0.001)
    boxplot(lapply(datRan, function(x) x[ , dd]), ylab = NA, axes = FALSE, lwd = 2, ylim = rr,
            bboxlwd = 6, boxwex = 0.4, whisklty = 1, whisklwd = 3,
            staplelty = 1, staplelwd = 3, medlwd = 4)
    yy <- axTicks(2)
    grid(lwd = 2)
    box(col = "grey40", lwd = 2)
    axis(side = 1, tck = -.01, labels = NA)
    mtext(c("MCMC", "VB", "CMC", "SDP", "WASP", "PIE"), at = 1:6, side = 1, line = 1, cex = 2)
    axis(side = 2, tck = -0.01, labels = NA, lwd = 1.5)
    mtext(yy, at = yy, side = 2, line = 0.3, las = 1, cex = 2)
    if(all(!(dd %in% c(4, 5)))) {
        mtext(bquote(Sigma[.(rnames[dd])]), side = 3, line = 0, adj = 0.5, cex = 2.5)
    } else {
        mtext(bquote(Sigma[.(rnames[dd])]), side = 3, line = 0, adj = 0.5, cex = 2.5)
    }
}
dev.off()


datCorr <- list(MCMC = mcmcCorr[[idxrf]][ , idxran],
               VB = vbCorr[[idxrf]][ , idxran],
               CMC = scotCorr[[idxrf]][ , idxran],
               SC = xingCorr[[idxrf]][ , idxran],
               WASP = waspCorr[[idxrf]][ , idxran],
               PIE = pieCorrEst[[idxrf]][ , idxran]
               )

pdf("~/pie/abrevaya/result/img/abe_ran_corr_box.pdf", 30, 10)
par(mfrow = c(2, 3))
par(cex = 1)
par(mar = c(2, 4, 3, 1), oma = c(1, 1, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (dd in 1:length(idxran)) {
    rr <- range(unlist(lapply(datCorr, function(x) x[ , dd]))) + c(-0.001, 0.001)
    boxplot(lapply(datCorr, function(x) x[ , dd]), ylab = NA, axes = FALSE, lwd = 2, ylim = rr,
            bboxlwd = 6, boxwex = 0.4, whisklty = 1, whisklwd = 3,
            staplelty = 1, staplelwd = 3, medlwd = 4)
    yy <- axTicks(2)
    grid(lwd = 2)
    box(col = "grey40", lwd = 2)
    axis(side = 1, tck = -.01, labels = NA)
    mtext(c("MCMC", "VB", "CMC", "SDP", "WASP", "PIE"), at = 1:6, side = 1, line = 1, cex = 2.5)
    if (dd == 6) {
    axis(side = 2, tck = -0.01, labels = NA, lwd = 1.5)
    mtext(round(yy, 1), at = round(yy, 1), side = 2, line = 0.3, las = 1, cex = 2)
    } else {
    axis(side = 2, tck = -0.01, labels = NA, lwd = 1.5)
    mtext(yy, at = yy, side = 2, line = 0.3, las = 1, cex = 2)
    }
    if (dd <=3) {
        mtext(bquote(Sigma[.(rnames[dd])]), side = 3, line = 0, adj = 0.5, cex = 2.5)
    } else {
        if (dd == 5) {
            mtext(bquote(rho[.(rnames[dd])]), side = 3, line = 0, adj = 0.5, cex = 2.5)
        } else {
            mtext(bquote(rho[.(rnames[dd])]), side = 3, line = 0, adj = 0.5, cex = 2.5)
        }
    }
}
dev.off()

fnames <- colnames(vbFix[[1]])
idxfix <- which(fnames %in% c("gestat", "male", "black", "adeqcode2", "adeqcode3", "pretri2"))

datFix <- list(MCMC = mcmcFix[[idxrf]][ , idxfix],
               MLE = mlFix[[idxrf]][ , idxfix],
               VB = vbFix[[idxrf]][ , idxfix],
               CMC = scotFix[[idxrf]][ , idxfix],
               SC = xingFix[[idxrf]][ , idxfix],
               WASP = waspFix[[idxrf]][ , idxfix],
               PIE = pieFixEst[[idxrf]][ , idxfix]
               )

pdf("~/pie/abrevaya/result/img/abe_fix_box_2row.pdf", 30, 10)
par(mfrow = c(2, 3))
par(cex = 1)
par(mar = c(2, 4, 2, 0), oma = c(1, 1, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (dd in 1:length(idxfix)) {
    rr <- range(unlist(lapply(datFix, function(x) x[ , dd]))) + c(-0.01, 0.01)
    boxplot(lapply(datFix, function(x) x[ , dd]), ylab = NA, axes = FALSE, lwd = 2, ylim = rr,
            bboxlwd = 6, boxwex = 0.4, whisklty = 1, whisklwd = 3,
            staplelty = 1, staplelwd = 3, medlwd = 4)
    yy <- axTicks(2)
    grid(lwd = 2)
    box(col = "grey40", lwd = 2)
    axis(side = 1, tck = -.01, labels = NA)
    mtext(c("MCMC", "  MLE", "VB", "CMC", "SDP", "WASP", "PIE"), at = 1:7, side = 1, line = 1, cex = 2)
    axis(side = 2, tck = -0.01, labels = NA, lwd = 1.5)
    mtext(round(yy, 2), at = yy, side = 2, line = 0.3, las = 1, cex = 2)
    mtext(fnames[idxfix[dd]], side = 3, line = 0, adj = 0.5, cex = 2.5)
}
dev.off()

resRan <- vector("list", 6)
names(resRan) <- c("MCMC", "VB", "PIE", "WASP", "Scot", "Xing")
for (jj in 1:6) {
    resRan[["MCMC"]][[jj]] <- list()
    resRan[["VB"]][[jj]] <- list()
    resRan[["PIE"]][[jj]] <- list()
    resRan[["WASP"]][[jj]] <- list()
    resRan[["Xing"]][[jj]] <- list()
    resRan[["Scot"]][[jj]] <- list()
    for (ii in 1:10) {
        rr <- range(c(mcmcRan[[ii]][ , jj], vbRan[[ii]][ , jj], pieRanEst[[ii]][ , jj],
                      waspRan[[ii]][ , jj]))
        bw1 <- dpik(mcmcRan[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw2 <- dpik(vbRan[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw3 <- dpik(pieRanEst[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw4 <- dpik(waspRan[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw5 <- dpik(scotRan[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw6 <- dpik(xingRan[[ii]][ , jj], range.x = rr, gridsize = 2000)
        dens1 <- bkde(mcmcRan[[ii]][ , jj], bandwidth = bw1, range.x = rr, gridsize = 2000)
        dens2 <- bkde(vbRan[[ii]][ , jj], bandwidth = bw2, range.x = rr, gridsize = 2000)
        dens3 <- bkde(pieRanEst[[ii]][ , jj], bandwidth = bw3, range.x = rr, gridsize = 2000)
        dens4 <- bkde(waspRan[[ii]][ , jj], bandwidth = bw4, range.x = rr, gridsize = 2000)
        dens5 <- bkde(scotRan[[ii]][ , jj], bandwidth = bw5, range.x = rr, gridsize = 2000)
        dens6 <- bkde(xingRan[[ii]][ , jj], bandwidth = bw6, range.x = rr, gridsize = 2000)
        resRan[["MCMC"]][[jj]][[ii]] <- dens1
        resRan[["VB"]][[jj]][[ii]] <- dens2
        resRan[["PIE"]][[jj]][[ii]] <- dens3
        resRan[["WASP"]][[jj]][[ii]] <- dens4
        resRan[["Scot"]][[jj]][[ii]] <- dens5
        resRan[["Xing"]][[jj]][[ii]] <- dens6
    }
}

resCorr <- vector("list", 6)
names(resCorr) <- c("MCMC", "VB", "PIE", "WASP", "Scot", "Xing")
for (jj in 1:6) {
    resCorr[["MCMC"]][[jj]] <- list()
    resCorr[["VB"]][[jj]] <- list()
    resCorr[["PIE"]][[jj]] <- list()
    resCorr[["WASP"]][[jj]] <- list()
    resCorr[["Xing"]][[jj]] <- list()
    resCorr[["Scot"]][[jj]] <- list()
    for (ii in 1:10) {
        rr <- range(c(mcmcCorr[[ii]][ , jj], vbCorr[[ii]][ , jj], pieCorrEst[[ii]][ , jj],
                      waspCorr[[ii]][ , jj]))
        bw1 <- dpik(mcmcCorr[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw2 <- dpik(vbCorr[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw3 <- dpik(pieCorrEst[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw4 <- dpik(waspCorr[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw5 <- dpik(scotCorr[[ii]][ , jj], range.x = rr, gridsize = 2000)
        bw6 <- dpik(xingCorr[[ii]][ , jj], range.x = rr, gridsize = 2000)
        dens1 <- bkde(mcmcCorr[[ii]][ , jj], bandwidth = bw1, range.x = rr, gridsize = 2000)
        dens2 <- bkde(vbCorr[[ii]][ , jj], bandwidth = bw2, range.x = rr, gridsize = 2000)
        dens3 <- bkde(pieCorrEst[[ii]][ , jj], bandwidth = bw3, range.x = rr, gridsize = 2000)
        dens4 <- bkde(waspCorr[[ii]][ , jj], bandwidth = bw4, range.x = rr, gridsize = 2000)
        dens5 <- bkde(scotCorr[[ii]][ , jj], bandwidth = bw5, range.x = rr, gridsize = 2000)
        dens6 <- bkde(xingCorr[[ii]][ , jj], bandwidth = bw6, range.x = rr, gridsize = 2000)
        resCorr[["MCMC"]][[jj]][[ii]] <- dens1
        resCorr[["VB"]][[jj]][[ii]] <- dens2
        resCorr[["PIE"]][[jj]][[ii]] <- dens3
        resCorr[["WASP"]][[jj]][[ii]] <- dens4
        resCorr[["Scot"]][[jj]][[ii]] <- dens5
        resCorr[["Xing"]][[jj]][[ii]] <- dens6
    }
}

resFix <- vector("list", 7)
names(resFix) <- c("MCMC", "ML", "VB", "PIE", "WASP", "Scot", "Xing")
for (jj in 1:14) {
    resFix[["MCMC"]][[jj]] <- list()
    resFix[["ML"]][[jj]] <- list()
    resFix[["VB"]][[jj]] <- list()
    resFix[["PIE"]][[jj]] <- list()
    resFix[["WASP"]][[jj]] <- list()
    resFix[["Scot"]][[jj]] <- list()
    resFix[["Xing"]][[jj]] <- list()
    for (ii in 1:10) {
        rr <- range(c(mcmcFix[[ii]][ , jj], vbFix[[ii]][ , jj], pieFixEst[[ii]][ , jj], waspFix[[ii]][ , jj],
                      mlFix[[ii]][ , jj], scotFix[[ii]][ , jj], xingFix[[ii]][ , jj]))
        bw1 <- dpik(mcmcFix[[ii]][ , jj], range.x = rr)
        bw2 <- dpik(vbFix[[ii]][ , jj], range.x = rr)
        bw3 <- dpik(pieFixEst[[ii]][ , jj], range.x = rr)
        bw4 <- dpik(waspFix[[ii]][ , jj], range.x = rr)
        bw5 <- dpik(mlFix[[ii]][ , jj], range.x = rr)
        bw6 <- dpik(scotFix[[ii]][ , jj], range.x = rr)
        bw7 <- dpik(xingFix[[ii]][ , jj], range.x = rr)
        dens1 <- bkde(mcmcFix[[ii]][ , jj], bandwidth = bw1, range.x = rr)
        dens2 <- bkde(vbFix[[ii]][ , jj], bandwidth = bw2, range.x = rr)
        dens3 <- bkde(pieFixEst[[ii]][ , jj], bandwidth = bw3, range.x = rr)
        dens4 <- bkde(waspFix[[ii]][ , jj], bandwidth = bw4, range.x = rr)
        dens5 <- bkde(mlFix[[ii]][ , jj], bandwidth = bw5, range.x = rr)
        dens6 <- bkde(scotFix[[ii]][ , jj], bandwidth = bw6, range.x = rr)
        dens7 <- bkde(xingFix[[ii]][ , jj], bandwidth = bw7, range.x = rr)
        resFix[["MCMC"]][[jj]][[ii]] <- dens1
        resFix[["VB"]][[jj]][[ii]] <- dens2
        resFix[["PIE"]][[jj]][[ii]] <- dens3
        resFix[["WASP"]][[jj]][[ii]] <- dens4
        resFix[["ML"]][[jj]][[ii]] <- dens5
        resFix[["Scot"]][[jj]][[ii]] <- dens6
        resFix[["Xing"]][[jj]][[ii]] <- dens7
    }
}

accRan <- array(0.0, dim = c(10, 6, 5))
for (jj in 1:6) {
    for (ii in 1:10) {
        accRan[ii, jj, 1] <- max(1 - sum(abs(resRan[["MCMC"]][[jj]][[ii]]$y  - resRan[["VB"]][[jj]][[ii]]$y) * diff(resRan[["VB"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accRan[ii, jj, 2] <- max(1 - sum(abs(resRan[["MCMC"]][[jj]][[ii]]$y - resRan[["PIE"]][[jj]][[ii]]$y) * diff(resRan[["PIE"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accRan[ii, jj, 3] <- max(1 - sum(abs(resRan[["MCMC"]][[jj]][[ii]]$y - resRan[["WASP"]][[jj]][[ii]]$y) * diff(resRan[["WASP"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accRan[ii, jj, 4] <- max(1 - sum(abs(resRan[["MCMC"]][[jj]][[ii]]$y - resRan[["Scot"]][[jj]][[ii]]$y) * diff(resRan[["Scot"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accRan[ii, jj, 5] <- max(1 - sum(abs(resRan[["MCMC"]][[jj]][[ii]]$y - resRan[["Xing"]][[jj]][[ii]]$y) * diff(resRan[["Xing"]][[jj]][[ii]]$x)[1]) / 2, 0)
    }
}

accCorr <- array(0.0, dim = c(10, 6, 5))
for (jj in 1:6) {
    for (ii in 1:10) {
        accCorr[ii, jj, 1] <- max(1 - sum(abs(resCorr[["MCMC"]][[jj]][[ii]]$y  - resCorr[["VB"]][[jj]][[ii]]$y) * diff(resCorr[["VB"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accCorr[ii, jj, 2] <- max(1 - sum(abs(resCorr[["MCMC"]][[jj]][[ii]]$y - resCorr[["PIE"]][[jj]][[ii]]$y) * diff(resCorr[["PIE"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accCorr[ii, jj, 3] <- max(1 - sum(abs(resCorr[["MCMC"]][[jj]][[ii]]$y - resCorr[["WASP"]][[jj]][[ii]]$y) * diff(resCorr[["WASP"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accCorr[ii, jj, 4] <- max(1 - sum(abs(resCorr[["MCMC"]][[jj]][[ii]]$y - resCorr[["Scot"]][[jj]][[ii]]$y) * diff(resCorr[["Scot"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accCorr[ii, jj, 5] <- max(1 - sum(abs(resCorr[["MCMC"]][[jj]][[ii]]$y - resCorr[["Xing"]][[jj]][[ii]]$y) * diff(resCorr[["Xing"]][[jj]][[ii]]$x)[1]) / 2, 0)
    }
}

accFix <- array(0.0, dim = c(10, 14, 6))
for (jj in 1:14) {
    for (ii in 1:10) {
        accFix[ii, jj, 1] <- max(1 - sum(abs(resFix[["MCMC"]][[jj]][[ii]]$y - resFix[["ML"]][[jj]][[ii]]$y) * diff(resFix[["ML"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accFix[ii, jj, 2] <- max(1 - sum(abs(resFix[["MCMC"]][[jj]][[ii]]$y - resFix[["VB"]][[jj]][[ii]]$y) * diff(resFix[["VB"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accFix[ii, jj, 3] <- max(1 - sum(abs(resFix[["MCMC"]][[jj]][[ii]]$y - resFix[["PIE"]][[jj]][[ii]]$y) * diff(resFix[["PIE"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accFix[ii, jj, 4] <- max(1 - sum(abs(resFix[["MCMC"]][[jj]][[ii]]$y - resFix[["WASP"]][[jj]][[ii]]$y) * diff(resFix[["WASP"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accFix[ii, jj, 5] <- max(1 - sum(abs(resFix[["MCMC"]][[jj]][[ii]]$y - resFix[["Scot"]][[jj]][[ii]]$y) * diff(resFix[["Scot"]][[jj]][[ii]]$x)[1]) / 2, 0)
        accFix[ii, jj, 6] <- max(1 - sum(abs(resFix[["MCMC"]][[jj]][[ii]]$y - resFix[["Xing"]][[jj]][[ii]]$y) * diff(resFix[["Xing"]][[jj]][[ii]]$x)[1]) / 2, 0)
    }
}

t(round(apply(accFix, 3, colMeans), 2))
t(round(apply(accFix, 3, colSds), 2))
tblFix <- matrix(paste(format(t(round(apply(accFix, 3, colMeans), 2)), nsmall = 2),
      paste0("(", format(t(round(apply(accFix, 3, colSds), 2)), nsmall = 2), ")")
      ), 6, 14)
rownames(tblFix) <- c("ML", "VB", "PIE", "WASP", "CMC", "SC")
xtable::xtable(tblFix)

tblRan <- matrix(paste(format(t(round(apply(accRan, 3, colMeans), 2)), nsmall = 2),
      paste0("(", format(t(round(apply(accRan, 3, colSds), 2)), nsmall = 2), ")")
      ), 5, 6)
rownames(tblRan) <- c("VB", "PIE", "WASP", "CMC", "SC")
xtable::xtable(tblRan)

tblCorr <- matrix(paste(format(t(round(apply(accCorr[ , c(1, 4, 6, 2, 3, 5), ], 3, colMeans), 2)), nsmall = 2),
      paste0("(", format(t(round(apply(accCorr[ , c(1, 4, 6, 2, 3, 5), ], 3, colSds), 2)), nsmall = 2), ")")
      ), 5, 6)
rownames(tblCorr) <- c("VB", "PIE", "WASP", "CMC", "SDP")
xtable::xtable(tblCorr)

fnames <- colnames(vbFix[[1]])

for (jj in 2:14) {
    pdf(paste0("~/pie/abrevaya/result/img/fix_", fnames[jj], ".pdf"), 40, 15)
    par(mfrow = c(2, 5))
    par(cex = 1)
    par(mar = c(0, 0, 0, 0), oma = c(4.5, 4.5, 0.4, 0.4))
    par(tcl = -0.02)
    par(mgp = c(2, 0.6, 0))
    for (ii in 1:10) {
        rr <- round(range(unlist(lapply(resFix[["MCMC"]][[jj]], function(x) x$y)),
                          unlist(lapply(resFix[["ML"]][[jj]], function(x) x$y)),
                          unlist(lapply(resFix[["VB"]][[jj]], function(x) x$y)),
                          unlist(lapply(resFix[["PIE"]][[jj]], function(x) x$y)),
                          unlist(lapply(resFix[["WASP"]][[jj]], function(x) x$y))))
        plot(resFix[["MCMC"]][[jj]][[ii]]$x, resFix[["MCMC"]][[jj]][[ii]]$y, col = colors[1], ylab = NA, xlab = NA, axes = FALSE, type = "l", ylim = rr)
        xxlab <- axTicks(1)
        yylab <- axTicks(2)
        if (ii == 1 | ii == 6) {
            axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
            mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 3, las = 2)
        }
        if (ii > 5) {
            axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
            mtext(format(round(xxlab, 2)), at = round(xxlab, 2), side = 1, line = 1, cex = 2.5)
        }
        grid(lwd = 5)
        box(col = "grey40", lwd = 5)
        lines(resFix[["MCMC"]][[jj]][[ii]]$x, resFix[["MCMC"]][[jj]][[ii]]$y, col = colors[1], lwd = 6)
        lines(resFix[["ML"]][[jj]][[ii]]$x, resFix[["ML"]][[jj]][[ii]]$y, col = colors[5], lwd = 6)
        lines(resFix[["VB"]][[jj]][[ii]]$x, resFix[["VB"]][[jj]][[ii]]$y, col = colors[2], lwd = 6)
        lines(resFix[["PIE"]][[jj]][[ii]]$x, resFix[["PIE"]][[jj]][[ii]]$y, col = colors[3], lwd = 6)
        lines(resFix[["WASP"]][[jj]][[ii]]$x, resFix[["WASP"]][[jj]][[ii]]$y, col = colors[4], lwd = 6)
        if (ii == 1) {
            legend("topright", c("MCMC", "ML", "VB", "PIE", "WASP"), col = colors[c(1, 5, 2:4)],
                   lty = "solid", lwd = 6, bty = "n", cex = 2)
        }
    }
    dev.off()
}

rnames <- c("dmage_dmage",  "nlbnl_nlbnl", "gestat_gestat")
cnames <- c("dmage_nbnl",  "dmage_gestat", "nlbnl_gestat")
varMap <- c(1, 4, 6)
covMap <- c(2, 3, 5)

for (jj in 1:3) {
    pdf(paste0("~/pie/abrevaya/result/img/var_", rnames[jj], ".pdf"), 40, 15)
    par(mfrow = c(2, 5))
    par(cex = 1)
    par(mar = c(0, 0, 0, 0), oma = c(4.5, 4.5, 0.4, 0.4))
    par(tcl = -0.02)
    par(mgp = c(2, 0.6, 0))
    for (ii in 1:10) {
        rr0 <- (range(unlist(lapply(resRan[["MCMC"]][[varMap[jj]]], function(x) x$y)), unlist(lapply(resRan[["PIE"]][[varMap[jj]]], function(x) x$y)), unlist(lapply(resRan[["WASP"]][[varMap[jj]]], function(x) x$y))))
        rr <- c(round(rr0[1]), ceiling(rr0[2]))
        plot(resRan[["MCMC"]][[varMap[jj]]][[ii]]$x, resRan[["MCMC"]][[varMap[jj]]][[ii]]$y, col = colors[1], ylab = NA, xlab = NA, axes = FALSE, type = "l", ylim = rr)
        xxlab <- axTicks(1)
        yylab <- axTicks(2)
        if (ii == 1 | ii == 6) {
            axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
            mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 3, las = 2)
        }
        if (ii > 5) {
            axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
            mtext(format(round(xxlab, 2)), at = round(xxlab, 2), side = 1, line = 1, cex = 2.5)
        }
        grid(lwd = 5)
        box(col = "grey40", lwd = 5)
        lines(resRan[["MCMC"]][[varMap[jj]]][[ii]]$x, resRan[["MCMC"]][[varMap[jj]]][[ii]]$y, col = colors[1], lwd = 6)
        lines(resRan[["VB"]][[varMap[jj]]][[ii]]$x, resRan[["VB"]][[varMap[jj]]][[ii]]$y, col = colors[2], lwd = 6)
        lines(resRan[["PIE"]][[varMap[jj]]][[ii]]$x, resRan[["PIE"]][[varMap[jj]]][[ii]]$y, col = colors[3], lwd = 6)
        lines(resRan[["WASP"]][[varMap[jj]]][[ii]]$x, resRan[["WASP"]][[varMap[jj]]][[ii]]$y, col = colors[4], lwd = 6)
        abline(v = mlRan[[ii]][[varMap[jj]]], col = colors[5], lwd = 6)
        if (ii == 1) {
            legend("topright", c("MCMC", "ML", "VB", "PIE", "WASP"), col = colors[c(1, 5, 2:4)],
                   lty = "solid", lwd = 6, bty = "n", cex = 2)
        }
    }
    dev.off()
}

for (jj in 1:3) {
    pdf(paste0("~/pie/abrevaya/result/img/cov_", cnames[jj], ".pdf"), 40, 15)
    par(mfrow = c(2, 5))
    par(cex = 1)
    par(mar = c(0, 0, 0, 0), oma = c(4.5, 4.5, 0.4, 0.4))
    par(tcl = -0.02)
    par(mgp = c(2, 0.6, 0))
    for (ii in 1:10) {
        rr0 <- (range(unlist(lapply(resRan[["MCMC"]][[covMap[jj]]], function(x) x$y)), unlist(lapply(resRan[["PIE"]][[covMap[jj]]], function(x) x$y)), unlist(lapply(resRan[["WASP"]][[covMap[jj]]], function(x) x$y))))
        rr <- c(round(rr0[1]), ceiling(rr0[2]))
        plot(resRan[["MCMC"]][[covMap[jj]]][[ii]]$x, resRan[["MCMC"]][[covMap[jj]]][[ii]]$y, col = colors[1], ylab = NA, xlab = NA, axes = FALSE, type = "l", ylim = rr)
        xxlab <- axTicks(1)
        yylab <- axTicks(2)
        if (ii == 1 | ii == 6) {
            axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
            mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 3, las = 2)
        }
        if (ii > 5) {
            axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
            mtext(format(round(xxlab, 2)), at = round(xxlab, 2), side = 1, line = 1, cex = 2.5)
        }
        grid(lwd = 5)
        box(col = "grey40", lwd = 5)
        lines(resRan[["MCMC"]][[covMap[jj]]][[ii]]$x, resRan[["MCMC"]][[covMap[jj]]][[ii]]$y, col = colors[1], lwd = 6)
        lines(resRan[["VB"]][[covMap[jj]]][[ii]]$x, resRan[["VB"]][[covMap[jj]]][[ii]]$y, col = colors[2], lwd = 6)
        lines(resRan[["PIE"]][[covMap[jj]]][[ii]]$x, resRan[["PIE"]][[covMap[jj]]][[ii]]$y, col = colors[3], lwd = 6)
        lines(resRan[["WASP"]][[covMap[jj]]][[ii]]$x, resRan[["WASP"]][[covMap[jj]]][[ii]]$y, col = colors[4], lwd = 6)
        abline(v = mlRan[[ii]][[covMap[jj]]], col = colors[5], lwd = 6)
        if (ii == 1) {
            legend("topright", c("MCMC", "ML", "VB", "PIE", "WASP"), col = colors[c(1, 5, 2:4)],
                   lty = "solid", lwd = 6, bty = "n", cex = 2)
        }
    }
    dev.off()
}


runTime <- list(log10(unlist(mcmcTime)),
                log10(unlist(vbTime)),
                log10(unlist(scotTime)),
                log10(unlist(xingTime)),
                log10(unlist(waspTime)),
                log10(unlist(pieTime))
                )

pdf("~/pie/abrevaya/result/img/abe_time.pdf", 8, 8)
par(cex = 1)
par(mar = c(0, 0, 0, 0), oma = c(3, 5.5, 0.2, 0.2))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
boxplot(runTime, ylab = NA, axes = FALSE, lwd = 2, ylim = c(1, 5.5))
grid(lwd=3)
box(col = "grey40", lwd = 4)
axis(side = 1, tck = -.01, labels = NA)
mtext(c("MCMC", "VB", "CMC", "SDP", "WASP", "PIE"), at = 1:6, side = 1, line = 1, cex = 2)
axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
mtext(format(seq(1, 5.5, by = 0.5)), at = seq(1, 5.5, by =0.5), side = 2, line = 0.3, las = 1, cex = 2)
mtext(expression(log[10] * " Seconds"), side = 2, outer = TRUE, cex = 2.5, line = 3)
dev.off()


ciRan <- array(0.0, dim = c(10, 6, 2, 6),
               dimnames = list(paste0("cv", 1:10),
                               paste0("dim", 1:6),
                               c("lb", "ub"),
                               c("MCMC", "VB", "PIE", "WASP", "Scot", "Xing")))
covMap <- c(1, 4, 6, 2, 3, 5)
for (ii in 1:10) {
    ciRan[ii, , , 1] <- do.call(rbind,
                                  lapply(split(mcmcRan[[ii]][ , covMap], col(mcmcRan[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciRan[ii, , , 2] <- do.call(rbind,
                                lapply(split(vbRan[[ii]][ , covMap], col(vbRan[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciRan[ii, , , 3] <- do.call(rbind,
                                lapply(split(pieRanEst[[ii]][ , covMap], col(pieRanEst[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciRan[ii, , , 4] <- do.call(rbind,
                                lapply(split(waspRan[[ii]][ , covMap], col(waspRan[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciRan[ii, , , 5] <- do.call(rbind,
                                lapply(split(scotRan[[ii]][ , covMap], col(scotRan[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciRan[ii, , , 6] <- do.call(rbind,
                                lapply(split(xingRan[[ii]][ , covMap], col(xingRan[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
}


tb2 <- apply(ciRan, c(2, 4), colMeans)
ciRanTbl <- matrix("0", 6, 6)
rownames(ciRanTbl) <- c("MCMC", "VB", "PIE", "WASP", "Scot", "Xing")
for (ii in 1:6) {
    for (jj in 1:6) {
        ciRanTbl[jj , ii] <- paste0("(",
                               format(round(tb2["lb", ii, jj], 4), nsmall = 4), ", ",
                               format(round(tb2["ub", ii, jj], 4), nsmall = 4), ")")
    }
}

xtable(ciRanTbl)

ciCorr <- array(0.0, dim = c(10, 6, 2, 6),
               dimnames = list(paste0("cv", 1:10),
                               paste0("dim", 1:6),
                               c("lb", "ub"),
                               c("MCMC", "VB", "PIE", "WASP", "Scot", "Xing")))
covMap <- c(1, 4, 6, 2, 3, 5)
for (ii in 1:10) {
    ciCorr[ii, , , 1] <- do.call(rbind,
                                  lapply(split(mcmcCorr[[ii]][ , covMap], col(mcmcCorr[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciCorr[ii, , , 2] <- do.call(rbind,
                                lapply(split(vbCorr[[ii]][ , covMap], col(vbCorr[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciCorr[ii, , , 3] <- do.call(rbind,
                                lapply(split(pieCorrEst[[ii]][ , covMap], col(pieCorrEst[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciCorr[ii, , , 4] <- do.call(rbind,
                                lapply(split(waspCorr[[ii]][ , covMap], col(waspCorr[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciCorr[ii, , , 5] <- do.call(rbind,
                                lapply(split(scotCorr[[ii]][ , covMap], col(scotCorr[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciCorr[ii, , , 6] <- do.call(rbind,
                                lapply(split(xingCorr[[ii]][ , covMap], col(xingCorr[[ii]][ , covMap])), function(x) quantile(x, prob = c(0.05, 0.95))))
}

tb2 <- apply(ciCorr, c(2, 4), colMeans)
ciCorrTbl <- matrix("0", 6, 6)
rownames(ciCorrTbl) <- c("MCMC", "VB", "PIE", "WASP", "Scot", "Xing")
for (ii in 1:6) {
    for (jj in 1:6) {
        ciCorrTbl[jj , ii] <- paste0("(",
                               format(round(tb2["lb", ii, jj], 3), nsmall = 3), ", ",
                               format(round(tb2["ub", ii, jj], 3), nsmall = 3), ")")
    }
}

xtable::xtable(ciCorrTbl)

ciFix <- array(0.0, dim = c(10, 14, 2, 7),
               dimnames = list(paste0("cv", 1:10),
                               paste0("dim", 1:14),
                               c("lb", "ub"),
                               c("MCMC", "ML", "VB", "PIE", "WASP", "Scot", "Xing")))

for (ii in 1:10) {
    ciFix[ii, , , 1] <- do.call(rbind,
                                lapply(split(mcmcFix[[ii]], col(mcmcFix[[ii]])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciFix[ii, , , 2] <- do.call(rbind,
                                lapply(split(mlFix[[ii]], col(mlFix[[ii]])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciFix[ii, , , 3] <- do.call(rbind,
                                lapply(split(vbFix[[ii]], col(vbFix[[ii]])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciFix[ii, , , 4] <- do.call(rbind,
                                lapply(split(pieFixEst[[ii]], col(pieFixEst[[ii]])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciFix[ii, , , 5] <- do.call(rbind,
                                lapply(split(waspFix[[ii]], col(waspFix[[ii]])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciFix[ii, , , 6] <- do.call(rbind,
                                lapply(split(scotFix[[ii]], col(scotFix[[ii]])), function(x) quantile(x, prob = c(0.05, 0.95))))
    ciFix[ii, , , 7] <- do.call(rbind,
                                lapply(split(xingFix[[ii]], col(xingFix[[ii]])), function(x) quantile(x, prob = c(0.05, 0.95))))
}

tb1 <- apply(ciFix, c(2, 4), colMeans)
ciFixTbl <- matrix("0", 7, 14)
for (ii in 1:14) {
    for (jj in 1:7) {
        ciFixTbl[jj , ii] <- paste0("(",
                               format(round(tb1["lb", ii, jj], 2), nsmall = 2), ", ",
                               format(round(tb1["ub", ii, jj], 2), nsmall = 2), ")")
    }
}

rownames(ciFixTbl) <- c("MCMC", "ML", "VB", "PIE", "WASP", "Scot", "Xing")

est <- list()
for(ii in 1:5) {
    est[[ii]] <- list()
}

for (ii in 1:10) {
    est[[1]][[ii]] <- do.call(rbind,
                              lapply(split(mcmcFix[[ii]], col(mcmcFix[[ii]])),
                                     function(x) quantile(x, prob = c(0.05, 0.95))))
    est[[2]][[ii]] <- do.call(rbind,
                              lapply(split(mlFix[[ii]], col(mlFix[[ii]])),
                                     function(x) quantile(x, prob = c(0.05, 0.95))))
    est[[3]][[ii]] <- do.call(rbind,
                              lapply(split(vbFix[[ii]], col(vbFix[[ii]])),
                                     function(x) quantile(x, prob = c(0.05, 0.95))))
    est[[4]][[ii]] <- do.call(rbind,
                              lapply(split(pieFixEst[[ii]], col(pieFixEst[[ii]])),
                                     function(x) quantile(x, prob = c(0.05, 0.95))))
    est[[5]][[ii]] <- do.call(rbind,
                              lapply(split(waspFix[[ii]], col(waspFix[[ii]])),
                                     function(x) quantile(x, prob = c(0.05, 0.95))))
}


xtable::xtable(res)

est <- list()
for(ii in 1:5) {
    est[[ii]] <- list()
}

for (ii in 1:10) {
    est[[1]][[ii]] <- do.call(rbind,
                              lapply(split(mcmcRan[[ii]], col(mcmcRan[[ii]])),
                                     function(x) quantile(x, prob = c(0.05, 0.95))))
    est[[2]][[ii]] <- cbind(mlRan[[ii]], mlRan[[ii]])
    est[[3]][[ii]] <- do.call(rbind,
                              lapply(split(vbRan[[ii]], col(vbRan[[ii]])),
                                     function(x) quantile(x, prob = c(0.05, 0.95))))
    est[[4]][[ii]] <- do.call(rbind,
                              lapply(split(pieRanEst[[ii]], col(pieRanEst[[ii]])),
                                     function(x) quantile(x, prob = c(0.05, 0.95))))
    est[[5]][[ii]] <- do.call(rbind,
                              lapply(split(waspRan[[ii]], col(waspRan[[ii]])),
                                     function(x) quantile(x, prob = c(0.05, 0.95))))
}

res <- matrix("0", 5, 6)
covMap <- c(1, 4, 6, 2, 3, 5)
for (ii in 1:6) {
    for (jj in 1:5) {
        res[jj , ii] <- paste0("(",
                               round(colMeans(do.call(rbind, lapply(est[[jj]], function(x) x[, 1])))[covMap[ii]], 2), ", ",
                               round(colMeans(do.call(rbind, lapply(est[[jj]], function(x) x[, 2])))[covMap[ii]], 2), ")")
    }
}

xtable::xtable(res)

## plots

fnames <- colnames(vbFix[[1]])

pdf(paste0("~/pie/abrevaya/result/img/abe_fixef_plot.pdf"), 60, 10)
par(mfrow = c(1, 6))
par(cex = 1)
par(mar = c(4, 4, 0, 0), oma = c(1.5, 1.5, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (jj in c(4, 5, 9, 10, 11, 13)) {
    for (ii in 5) {
        rr0 <- c(range(unlist(lapply(resFix[["MCMC"]][[jj]], function(x) x$y)),
                       unlist(lapply(resFix[["ML"]][[jj]], function(x) x$y)),
                       unlist(lapply(resFix[["VB"]][[jj]], function(x) x$y)),
                       unlist(lapply(resFix[["PIE"]][[jj]], function(x) x$y)),
                       unlist(lapply(resFix[["WASP"]][[jj]], function(x) x$y))))
        rr <- c(0, ceiling(rr0[2]) + 0.5)
        plot(resFix[["MCMC"]][[jj]][[ii]]$x, resFix[["MCMC"]][[jj]][[ii]]$y, col = colors[1], ylab = NA, xlab = NA, axes = FALSE, type = "l", ylim = rr)
        xxlab <- axTicks(1)
        yylab <- axTicks(2)
        axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 3, las = 2)
        axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(round(xxlab, 2)), at = round(xxlab, 2), side = 1, line = 1, cex = 3)
        grid(lwd = 5)
        box(col = "grey40", lwd = 5)
        mtext(fnames[jj], side = 3, line = -3.5, adj = 0.1, cex = 5)
        lines(resFix[["MCMC"]][[jj]][[ii]]$x, resFix[["MCMC"]][[jj]][[ii]]$y, col = colors[1], lwd = 6)
        lines(resFix[["ML"]][[jj]][[ii]]$x, resFix[["ML"]][[jj]][[ii]]$y, col = colors[5], lwd = 6)
        lines(resFix[["VB"]][[jj]][[ii]]$x, resFix[["VB"]][[jj]][[ii]]$y, col = colors[2], lwd = 6)
        lines(resFix[["WASP"]][[jj]][[ii]]$x, resFix[["WASP"]][[jj]][[ii]]$y, col = colors[4], lwd = 6)
        lines(resFix[["PIE"]][[jj]][[ii]]$x, resFix[["PIE"]][[jj]][[ii]]$y, col = colors[3], lwd = 6)
        if (jj == 4) {
            legend("topright", c("MCMC", "ML", "VB", "WASP", "PIE"), col = colors[c(1, 5, 2, 4, 3)],
                   lty = "solid", lwd = 6, bty = "n", cex = 2)
        }
    }
}
dev.off()


rnames <- c("(dmage, dmage)",  "(nlbnl, nlbnl)", "(gestat, gestat)",
            "(dmage, nbnl)",  "(dmage, gestat)", "(nlbnl, gestat)")
vmap <- c(1, 4, 6, 2, 3, 5)

pdf("~/pie/abrevaya/result/img/abe_ranef_plot.pdf", 60, 10)
par(mfrow = c(1, 6))
par(cex = 1)
par(mar = c(0, 6, 0, 0), oma = c(4.5, 0.4, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (jj in 1:6) {
    for (ii in 5) {
        rr0 <- (range(unlist(lapply(resRan[["MCMC"]][[vmap[jj]]], function(x) x$y)), unlist(lapply(resRan[["PIE"]][[vmap[jj]]], function(x) x$y)), unlist(lapply(resRan[["WASP"]][[vmap[jj]]], function(x) x$y))))
        rr <- c(round(rr0[1]), ceiling(rr0[2]))
        plot(resRan[["MCMC"]][[vmap[jj]]][[ii]]$x, resRan[["MCMC"]][[vmap[jj]]][[ii]]$y, col = colors[1], ylab = NA, xlab = NA, axes = FALSE, type = "l", ylim = rr)
        xxlab <- axTicks(1)
        yylab <- axTicks(2)
        axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 3, las = 2)
        axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(round(xxlab, 2)), at = round(xxlab, 2), side = 1, line = 1.2, cex = 3)
        grid(lwd = 5)
        box(col = "grey40", lwd = 5)
        mtext(rnames[jj], side = 3, line = -3.5, adj = 0.1, cex = 4)
        lines(resRan[["MCMC"]][[vmap[jj]]][[ii]]$x, resRan[["MCMC"]][[vmap[jj]]][[ii]]$y, col = colors[1], lwd = 6)
        lines(resRan[["VB"]][[vmap[jj]]][[ii]]$x, resRan[["VB"]][[vmap[jj]]][[ii]]$y, col = colors[2], lwd = 6)
        lines(resRan[["WASP"]][[vmap[jj]]][[ii]]$x, resRan[["WASP"]][[vmap[jj]]][[ii]]$y, col = colors[4], lwd = 6)
        lines(resRan[["PIE"]][[vmap[jj]]][[ii]]$x, resRan[["PIE"]][[vmap[jj]]][[ii]]$y, col = colors[3], lwd = 6)
        abline(v = mlRan[[ii]][[vmap[jj]]], col = colors[5], lwd = 6)
        if (jj == 1) {
            legend("topright", c("MCMC", "ML", "VB", "WASP", "PIE"), col = colors[c(1, 5, 2, 4, 3)],
                   lty = "solid", lwd = 6, bty = "n", cex = 2)
        }
    }
}
dev.off()
