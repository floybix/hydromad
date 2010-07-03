## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


fitBySCE <-
    function(MODEL, 
             objective = hydromad.getOption("objective"),
             control = hydromad.getOption("sce.control"),
             vcov = FALSE)
{
    start_time <- proc.time()
    objective <- buildCachedObjectiveFun(objective, MODEL)
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
    if (isTRUE(hydromad.getOption("trace")))
            control$trace <- 1
    lower <- sapply(parlist, min)
    upper <- sapply(parlist, max)
    initpars <- sapply(parlist, mean) ## TODO: allow sampling?
    bestModel <- MODEL
    bestFunVal <- -Inf
    do_sce <- function(pars) {
        thisMod <- update(MODEL, newpars = pars)
        if (!isValidModel(thisMod))
            return(NA)
        thisVal <- objFunVal(thisMod, objective = objective)
        if (isTRUE(thisVal > bestFunVal)) {
            bestModel <<- thisMod
            bestFunVal <<- thisVal
        }
        ## SCEoptim does minimisation, so:
        return(- thisVal)
    }
    ans <- SCEoptim(do_sce, initpars, lower = lower, upper = upper,
                    control = control)
    if (ans$convergence != 0) {
        if (!isTRUE(hydromad.getOption("quiet"))) {
            warning(ans$message)
        }
        bestModel$msg <- ans$message
    }
    bestModel$funevals <- ans$counts
    bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
    bestModel$objective <- objective
    if (vcov) {
        ## estimate covariance matrix from final population
        ## TODO
        #bestModel$cov.mat <-
        #ans$POP.ALL
        warning("vcov not yet implemented")
    }
    bestModel$fit.call <- match.call()
    bestModel$fit.result <- ans
    return(bestModel)
}
