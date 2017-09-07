sampleFromMixMdl <- function (yvec, xmat, zmat, group, niter, nburn, nthin, dirp, id) {
    library(inline)
    library(Rcpp)
    library(rstan)

    gg <- ordered(as.character(group), levels = sort(unique(group)))
    group <- as.integer(gg)

    simList = list(
        nobs = length(yvec),
        nfixef = ncol(xmat),
        nranef = ncol(zmat),
        ngroup = length(unique(group)),
        xmat = xmat,
        zmat = zmat,
        group = group,
        yvec = yvec)

    seeds <- (1:1000) * as.numeric(gsub(":", "", substr(Sys.time(), 12, 19)))

    stanCode <- readChar("full_lme.stan", file.info("full_lme.stan")$size)
    startTime <- proc.time()
    mdl <- stan(model_code = stanCode, data = simList, iter = niter, warmup = nburn, chains = 1, thin = nthin, seed = seeds[id])
    endTime <- proc.time()

    lst <- mdl@sim$samples[[1]]
    bs <- grep("fixef|covRanef", names(lst))
    sampdf <- do.call(cbind, lst[bs])

    write.table(sampdf[(nrow(sampdf) - (niter - nburn) / nthin + 1):nrow(sampdf), ], paste0(dirp, "result/mcmc/samp_", id, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")

    saveRDS(list(samples = sampdf[(nrow(sampdf) - (niter - nburn) / nthin + 1):nrow(sampdf), ], time = endTime - startTime), paste0(dirp, "result/mcmc/samp_time_", id, ".rds"))

    return("done")
}
