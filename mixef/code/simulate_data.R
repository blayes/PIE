rm(list=ls())

set.seed(12345)

setwd("~/pie/mixef/code")
library(matrixStats)

## based on Perry (2017). JRSS-B.
## http://ptrckprry.com/code/
genData <- function (ngroup, nobs, nfixef = 4, nranef = 3, dispersion, family) {
    ## family
    if (is.character(family)) {
        family <- get(family, mode="function", envir=parent.frame())()
    } else if (is.function(family)) {
        family <- family()
    } else if (is.null(family$family)) {
        stop("'family' not recognized")
    }

    if (family$family %in% c("binomial", "poisson") && dispersion != 1) {
        warning(sprintf("dispersion parameter is ignored for '%s' family", family$family))
    }

    fixef <- rep(c(-1, 1), length = nfixef)

    if (nranef == 3) {
        ranef.corr <- matrix(c(1, -0.4, 0.3,
                               -0.4, 1, 0.001,
                               0.3, 0.001, 1),
                             nranef, nranef)

        ranef.cov <- outer(sqrt(1:nranef), sqrt(1:nranef)) * ranef.corr
        ranef.cov.sqrt <- chol(ranef.cov)
    } else {
        stop(sprintf("currently no default for '%s' rand effs", nranef))
    }

    ## generate coefficients
    u <- matrix(rnorm(ngroup * nranef), ngroup, nranef)
    ranef <- u %*% ranef.cov.sqrt

    ## generate group
    suppressWarnings({ # ignore warning about using Walker's alias method
        group <- sample.int(ngroup, nobs, replace=TRUE)
    })

    ## generate feature  matrices with Pr(x[i,j] = +1) = P(x[i,j] = -1) = 1/2,
    x <- matrix(sample(c(-1, +1), nobs * nfixef, replace=TRUE), nobs, nfixef)
    z <- matrix(sample(c(-1, +1), nobs * nranef, replace=TRUE), nobs, nranef)

    ## compute linear predictors and generate observations
    eta <- drop(x %*% fixef) + rowSums(z * ranef[group,])
    mu <- family$linkinv(eta)

    if (family$family == "gaussian") {
        y <- rnorm(nobs, mean=mu, sd=sqrt(dispersion))
    } else if (family$family == "binomial") {
        y <- as.numeric(rbinom(nobs, 1, mu))
    } else if (family$family == "poisson") {
        y <- rpois(nobs, mu)
    } else {
        stop(sprintf("family '%s' not supported", family$family))
    }

    list(ngroup = ngroup, nobs = nobs,
         fixef = fixef, ranef = ranef,
         ranef.corr = ranef.corr,
         dispersion = dispersion,
         ranef.cov = ranef.cov,
         group = group, x = x, z = z, y.mean = mu, y = y, family=family)
}

ngroup <- 5000
nobs <- 100000

# fixed effects
nfixef <- 4

# random effects
nranef <- 3

# noise
dispersion <- 1

# family
family <- gaussian()

repData <- list()
for (ii in 1:10) {
    repData[[ii]] <- genData(ngroup = ngroup, nobs = nobs, nfixef = nfixef, nranef = nranef, dispersion = dispersion, family = family)
}

saveRDS(repData, paste0("../data/lme.rds"))

rm(list=ls())

set.seed(12345)

nrep <- 10
npart <- 20

reps <- readRDS("../data/lme.rds")
parts <- vector("list", length = nrep)
partsIdx <- vector("list", length = nrep)
for (ii in 1:nrep) {
    parts[[ii]] <- vector("list", length = npart)
}

for (r in seq_len(nrep)) {
  lst <- reps[[r]]
  grpSplit <- split(1:nrow(lst$x), lst$group)
  partsIdx <- sample(1:npart, length(grpSplit), replace = TRUE)
  for (ii in 1:npart) {
    grpIdx <- which(partsIdx == ii)
    idx <- unlist(grpSplit[grpIdx])
    parts[[r]][[ii]]$nobs <- length(idx)
    parts[[r]][[ii]]$x <- lst$x[idx, ]
    parts[[r]][[ii]]$y <- lst$y[idx]
    parts[[r]][[ii]]$z <- lst$z[idx, ]
    parts[[r]][[ii]]$group <- lst$group[idx]
    parts[[r]][[ii]]$y.mean <- lst$y.mean[idx]
    parts[[r]][[ii]]$ranef <- lst$ranef[lst$group[idx], ]
  }
}

saveRDS(parts, "../data/wasp_lme.rds")
