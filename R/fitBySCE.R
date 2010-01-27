## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2009 Felix Andrews <felix@nfrac.org>
##


fitBySCE <-
    function(MODEL, ...,
             objective = ihacres.getOption("objective"),
             control = list())
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
    ## remove any fixed parameters
    parlist <- parlist[!isfixed]
    if (isTRUE(ihacres.getOption("trace")))
            control$trace <- 1
    lower <- sapply(parlist, min)
    upper <- sapply(parlist, max)
    initpars <- sapply(parlist, mean) ## TODO: allow sampling?
    bestModel <- MODEL
    bestFunVal <- Inf
    do_sce <- function(pars) {
        thisMod <- do.call("update", c(quote(MODEL), as.list(pars)))
        if (!isValidModel(thisMod))
            return(NA)
        thisVal <- objFunVal(thisMod, objective = objective)
        if (thisVal < bestFunVal) {
            bestModel <<- thisMod
            bestFunVal <<- thisVal
        }
        thisVal
    }
    ans <- SCEoptim(do_sce, initpars, lower = lower, upper = upper,
                    control = control)
    if (ans$convergence != 0) {
        if (!isTRUE(ihacres.getOption("quiet"))) {
            warning(ans$message)
        }
        bestModel$msg <- ans$message
    }
    bestModel$funevals <- ans$counts
    bestModel$timing <- proc.time() - start_time
    return(bestModel)
}
