## meanList = list of vector of multivariate mean
## covList = list of cov matrices
## wts = rep(1 / K, K), where K is the # of subset posteriors
obtainGaussianBarycenter <- function(meanList, covList, wts = NULL) {
    ncomp <- length(covList)
    ndim <- nrow(covList[[1]])
    baryMean <- rowMeans(do.call(cbind, meanList))

    if (is.null(wts)) wts <- rep(1, ncomp) / ncomp

    # google: http://realizationsinbiostatistics.blogspot.com/2008/08/matrix-square-roots-in-r_18.html
    denman.beavers <- function(mat, maxit = 50) {
        stopifnot(nrow(mat) == ncol(mat))
        niter <- 0
        y <- mat
        z <- diag(rep(1,nrow(mat)))
        for (niter in 1:maxit) {
            y.temp <- 0.5*(y+solve(z))
            z <- 0.5*(z+solve(y))
            y <- y.temp
        }
        return(list(sqrt=y,sqrt.inv=z))
    }

    sqrtList <- lapply(covList, function(x) {
        sx <- denman.beavers(x)
        sx$sqrt
    })

    barySd <- sqrtList[[1]]
    for (ii in 2:ncomp) {
        barySd <- barySd + sqrtList[[ii]]
    }

    baryCov <- barySd %*% barySd
    err <- baryCov
    cnt <- 1
    while ((norm(err, type = "F") > 1e-10) & cnt < 100) {
        if (cnt %% 10 == 0)  cat("iter: ", cnt, "\n")
        sumMat <- matrix(0.0, nrow = ndim, ncol = ndim)
        for (ii in 1:ncomp) {
            tmp2 <- denman.beavers(barySd %*% covList[[ii]] %*% barySd)
            sumMat <- sumMat + wts[ii] * tmp2$sqrt
        }
        baryCov1 <- sumMat
        err <- baryCov - baryCov1
        baryCov <- baryCov1
        barySd <- denman.beavers(baryCov)$sqrt
        cnt <- cnt + 1
    }

    list(mean = baryMean, cov = baryCov)
}
