sampleGdp <- function (yvec, xmat, ngrid = 100, niter, nburn, nthin) {
  library(mgcv)
  library(MCMCpack)
  library(matrixStats)

  ndim <- ncol(xmat)
  nobs <- nrow(xmat)

  cts <- 0
  sigmaSamp <- rep(0.0, (niter - nburn) / nthin)
  alphaSamp <- rep(0.0, (niter - nburn) / nthin)
  etaSamp <- rep(0.0, (niter - nburn) / nthin)
  betaSamp <- matrix(0.0, (niter - nburn) / nthin, ndim)

  eta <- 1
  alpha <- 1
  betas <- runif(ndim, 1, 10)
  sigma2 <- 1
  lambdas <- runif(ndim, 1, 10)
  taus <- runif(ndim, 1, 10)
  agrid <- seq(0.01, 0.99, length = ngrid)
  logWtAlpha <- rep(0, length(agrid))
  nuErr <- 2; aErr <- 1; aA <- 10^2

  gramMat <- crossprod(xmat, xmat)

  startTime <- proc.time()
  for (its in 1:niter) {
      covBeta <- sigma2 * chol2inv(chol(gramMat + diag(1 / taus)))
      muBeta <- covBeta %*% crossprod(xmat, yvec) / sigma2
      betas <- as.numeric(muBeta + crossprod(chol(covBeta), rnorm(ndim)))

      ## update post. for hyper parameters in half t
      aErr <- rinvgamma(1, shape = 0.5 * (nuErr + 1), scale = nuErr / sigma2 + 1 / aA^2)
      ## update post. for error var.
      resids <- yvec - as.numeric(xmat %*% betas)
      shp <- 0.5 * (nobs + nuErr)
      scl <- sum(resids^2) / 2 + nuErr / aErr
      sigma2 <- rinvgamma(1, shape = shp, scale = scl)

      for (pp in 1:ndim) {
          lambdas[pp] <- rgamma(1, shape = alpha + 1, scale = eta + abs(betas[pp]) / sqrt(sigma2))
          taus[pp] <- 1 / rig(1, mean = abs(lambdas[pp] * sqrt(sigma2) / betas[pp]), scale = 1 / lambdas[pp]^2)
      }

      logWtAlpha <- ndim * log((1 - agrid) / agrid) - sum(log(1 + abs(betas) / (eta * sqrt(sigma2)))) / agrid
      probsAlpha <- exp(logWtAlpha - max(logWtAlpha))
      idx <- sample(seq_along(agrid), 1, prob = probsAlpha)
      alpha <- 1 / agrid[idx] - 1

      if (its %% 1000 == 0) cat("iteration: ", its, "\n")

      if (its > nburn & its %% nthin == 0) {
          cts <- cts + 1
          sigmaSamp[cts] <- sigma2
          betaSamp[cts, ] <- betas
          alphaSamp[cts] <- alpha
      }
  }
  endTime <- proc.time()

  list(
    sigmaSamp = sigmaSamp,
    betaSamp = betaSamp,
    alphaSamp = alphaSamp,
    time = endTime - startTime
    )
}
