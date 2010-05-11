## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


fitByDream <-
    function(MODEL,
             #func.type = "calc.rmse",
             #measurement = list(data = observed(MODEL)),
             loglik = ~ -0.5 * sum((Q-X)^2),
             control = hydromad.getOption("dream.control"),
             vcov = TRUE)
{
    library(dream)
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
    if (!isTRUE(hydromad.getOption("trace")))
            control$REPORT <- 0
    do_dream <- function(pars) {
        names(pars) <- names(parlist)
        thisMod <- update(MODEL, newpars = pars)
        if (!isValidModel(thisMod))
            return(-1e8)
        #as.numeric(fitted(thisMod))
        objFunVal(thisMod, objective = loglik)
    }
    ans <- dream(do_dream, pars = parlist,
                 func.type = "logposterior.density",
                 control = control)
    bestPars <- coef(ans, method = "sample.ml")
    bestModel <- update(MODEL, newpars = bestPars)
    bestModel$funevals <- ans$fun.evals
    bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
    bestModel$objective <- loglik
    if (vcov) {
        ## estimate covariance matrix from final population
        start <- end(ans$Sequences)/2+1
        bestModel$cov.mat <-
            cov(as.matrix(window(ans$Sequences, start = start)))
    }
    bestModel$fit.call <- match.call()
    bestModel$fit.result <- ans
    return(bestModel)
}
