## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2009 Felix Andrews <felix@nfrac.org>
##


fitBySampling <-
    function(MODEL,
             objective = ihacres.getOption("objective"),
             samples = ihacres.getOption("fit.samples"),
             sampletype = c("latin.hypercube", "random", "all.combinations"))
{
    start_time <- proc.time()
    parlist <- as.list(coef(MODEL, warn = FALSE))
    ## remove any missing parameters
    isok <- sapply(parlist, function(x) !any(is.na(x)))
    parlist <- parlist[isok]
    ## check which parameters are uniquely specified
    isfixed <- (sapply(parlist, length) == 1)
    if (all(isfixed)) {
        warning("all parameters are fixed, so can not fit")
        return(MODEL)
    }
    ## generate parameter sets
    psets <- parameterSets(parlist, samples = samples, method = sampletype)
    bestModel <- MODEL
    bestFunVal <- Inf
    for (i in seq(NROW(psets))) {
        if (isTRUE(ihacres.getOption("trace"))) {
            run_name <- paste(names(psets), format(psets[i,], digits=3),
                              sep = "=", collapse = ",")
            message(run_name)
        }
        thisPars <- as.list(psets[i,])
        thisMod <- do.call("update", c(quote(MODEL), thisPars))
        if (!isValidModel(thisMod))
            next
        thisVal <- objFunVal(thisMod, objective = objective)
        if (thisVal < bestFunVal) {
            bestModel <- thisMod
            bestFunVal <- thisVal
        }
    }
    bestModel$funevals <- NROW(psets)
    bestModel$timing <- proc.time() - start_time
    bestModel
}
