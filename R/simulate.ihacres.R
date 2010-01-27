## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


simulate.ihacres <-
    function(object, nsim, seed, ...,
             sampletype = c("random", "latin.hypercube", "all.combinations"),
             FUN = NULL)
{
    MODEL <- object
    samples <- nsim
    if (!missing(seed))
        set.seed(seed)
    if (is.character(FUN))
        FUN <- get(FUN, mode = "function")
    ## this is mostly a copy of fitBySampling
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
    result <- list()
    for (i in seq(NROW(psets))) {
        run_name <- paste(names(psets), format(psets[i,], digits=3),
                          sep = "=", collapse = ",")
        if (isTRUE(ihacres.getOption("trace"))) {
            message(run_name)
        }
        thisPars <- as.list(psets[i,])
        thisMod <- do.call("update", c(quote(MODEL), thisPars))
        ## store model or derived result
        result[[i]] <-
            if (is.null(FUN)) thisMod else FUN(thisMod, ...)
        names(result)[i] <- run_name
    }
    if (is.null(FUN))
        result <- as.runlist(result)
    attr(result, "psets") <- psets
    result
}

