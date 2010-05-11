## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


fitByOptim <-
    function(MODEL,
             objective = hydromad.getOption("objective"),
             method = hydromad.getOption("optim.method"),
             control = hydromad.getOption("optim.control"),
             samples = hydromad.getOption("fit.samples"),
             sampletype = c("latin.hypercube", "random", "all.combinations"),
             multistart = FALSE,
             vcov = FALSE,
             hessian = vcov,
             initpars = NULL)
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
    if (multistart) {
        ## multiple optimisation runs with different starting points.
        ## generate parameter sets
        psets <- parameterSets(parlist, samples = samples, method = sampletype)
        bestModel <- MODEL
        bestFunVal <- Inf
        objseq <- numeric()
        funevals <- 0
        for (i in seq(NROW(psets))) {
            thisPars <- as.list(psets[i,,drop=FALSE])
            thisMod <-
                try(fitByOptim(MODEL, objective = objective,
                               method = method, control = control,
                               hessian = hessian,
                               multistart = FALSE, initpars = thisPars),
                    silent = TRUE)
            if (!isValidModel(thisMod))
                next
            funevals <- funevals + thisMod$funevals
            objseq <- c(objseq, thisMod$fit.result$objseq)
            thisVal <- objFunVal(thisMod, objective = objective)
            if (thisVal < bestFunVal) {
                bestModel <- thisMod
                bestFunVal <- thisVal
            }
        }
        bestModel$funevals <- funevals
        bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
        bestModel$fit.call <- match.call()
        bestModel$fit.result$objseq <- objseq
        return(bestModel)
        
    } else {
        ## single optimisation run.
        lower <- sapply(parlist, min)
        upper <- sapply(parlist, max)
#        lower[is.na(lower)] <- -Inf
#        upper[is.na(upper)] <- Inf
        pre.funevals <- 0
        if (!is.null(initpars)) {
            ## initial parameter values were specified.
            initpars <- unlist(initpars)[names(parlist)]
        } else if (samples > 1) {
            ## do sampling to find a good starting point
            preMODEL <- fitBySampling(MODEL, objective = objective,
                                      samples = samples, sampletype = sampletype)
            if (!isValidModel(preMODEL))
                return(preMODEL)
            pre.funevals <- preMODEL$funevals
            initpars <- coef(preMODEL)[names(parlist)]
        } else {
            initpars <- sapply(parlist, mean)
        }
        ## now optimise
        control <- modifyList(hydromad.getOption("optim.control"),
                              control)
        if (isTRUE(hydromad.getOption("quiet")))
            control$trace <- 0
        bestModel <- MODEL
        bestFunVal <- Inf
        objseq <- rep(NA_real_, 100)
        i <- 0
        do_optim <- function(pars) {
            ## TODO: handle boundaries better
            #i <- which(pars < lower)
            #if (any(i)) {
            #    return( 1e12 + (sum(lower[i] - pars[i])) * 1e6 )
            #}
            #i <- which(pars > upper)
            #if (any(i)) {
            #    return( 1e12 + (sum(pars[i] - upper[i])) * 1e6 )
            #}

            i <<- i + 1
            ## TODO: just set params to bounded values?

            if (any(pars < lower)) return(NA)
            if (any(pars > upper)) return(NA)
            thisMod <- update(MODEL, newpars = pars)
            if (!isValidModel(thisMod))
                return(NA)
            thisVal <- objFunVal(thisMod, objective = objective)
            objseq[i] <<- thisVal
            if (thisVal < bestFunVal) {
                bestModel <<- thisMod
                bestFunVal <<- thisVal
            }
            thisVal
        }
        if (!isTRUE(hydromad.getOption("catch.errors.optim")))
            try <- force ## i.e. skip the try()
        lowerb <- if (method == "L-BFGS-B") lower else -Inf
        upperb <- if (method == "L-BFGS-B") upper else Inf
        ans <- try(optim(initpars, do_optim, method = method,
                         lower = lowerb, upper = upperb,
                         control = control, hessian = hessian))
        if (inherits(ans, "try-error")) {
            bestModel$msg <- ans
            return(bestModel)
        }
        if (ans$convergence != 0) {
            msg <- if (ans$convergence == 1) {
                "optim() reached maximum iterations"
            } else {
                paste("optim() returned convergence code",
                      ans$convergence)
            }
            if (!isTRUE(hydromad.getOption("quiet"))) {
                warning(msg)
            }
            bestModel$msg <- msg
        }
        bestModel$funevals <- ans$counts[1] + pre.funevals
        bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
        bestModel$objective <- objective
        if (vcov) {
            ## approximate covariance matrix from inverse of hessian (often poor!)
            bestModel$cov.mat <- solve(ans$hessian)
        }
        bestModel$fit.call <- match.call()
        ans$objseq <- objseq[1:i]
        bestModel$fit.result <- ans
        return(bestModel)
    }
}
