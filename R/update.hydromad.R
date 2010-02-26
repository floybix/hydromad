## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

absorbScale <- function(object, gain)
    UseMethod("absorbScale")

absorbScale.hydromad <- function(object, gain)
    return(NULL)

update.hydromad <-
    function(object, ..., newdata = NULL,
             sma, routing, rfit, warmup, weights,
             and.rescale = TRUE)
{
    ## update timestamp
    object$last.updated <- Sys.time()

    ## update the stored call
    upcall <- match.call()
    nm <- names(upcall)
    if (!is.null(newdata)) {
        object$call[2] <- upcall["newdata"]
    }
    if (!is.null(nm)) {
        nm <- nm[!(nm %in% c("", "object", "newdata", "and.rescale"))]
        if (length(nm) > 0) {
            extras <- as.list(upcall)[nm]
            ## round off parameter values in call (for display purposes)
            isnum <- sapply(extras, is.numeric)
            extras[isnum] <- lapply(extras[isnum], signif, 6)
            object$call[nm] <- extras
        }
    }
    ## update DATA
    RUNSMA <- is.null(object$U)
    if (!is.null(newdata)) {
        object$data <- as.ts(newdata)
        RUNSMA <- TRUE
    }
    ## update named items
    if (!missing(sma)) {
        if (!is.null(object$sma)) {
            ## remove old SMA-specific parameters
            object$parlist <- as.list(coef(object, "routing", warn = FALSE))
        }
        class(object) <- "hydromad"
        if (!is.null(sma)) {
            ## take default parameter ranges/values from hydromad.options()
            pardefs <- hydromad.getOption(sma)
            if (is.null(pardefs)) {
                pardef.fn <- paste(sma, ".ranges", sep = "")
                if (exists(pardef.fn, mode = "function"))
                    pardefs <- do.call(pardef.fn, list())
            }
            object$parlist <-
                modifyList(as.list(pardefs),
                           object$parlist)
            class(object) <- c(paste("hydromad", sma, sep="."),
                               "hydromad")
        }
        object$sma <- sma
        object$sma.fun <- NULL
        object$sma.args <- NULL
        if (!is.null(sma)) {
            object$sma.fun <- paste(sma, ".sim", sep = "")
            force(get(object$sma.fun, mode = "function")) ## check exists
            object$sma.args <- formals(object$sma.fun)
        }
        RUNSMA <- TRUE
    }
    if (!missing(routing)) {
        if (!is.null(object$routing)) {
            ## remove old routing-specific parameters
            object$parlist <- as.list(coef(object, "sma", warn = FALSE))
        }
        if (!is.null(routing)) {
            ## take default parameter ranges/values from hydromad.options()
            object$parlist <-
                modifyList(as.list(hydromad.getOption(routing)),
                           object$parlist)
        }
        object$routing <- routing
    }
    if (!missing(rfit))
        object$rfit <- rfit
    if (!missing(warmup)) {
        stopifnot(is.numeric(warmup))
        object$warmup <- warmup
    }
    if (!missing(weights)) {
        ## weights may be either a vector correpsonding to observed data,
        ## or a function to apply to observed data to calculate weights
        ## TODO: maybe better to allow a column "weights" in DATA?
        if (!is.null(weights) && !is.function(weights))
            stopifnot(NROW(weights) == NROW(object$data))
        object$weights <- weights
    }
    ## delete any previous stored results -- we will run again
    object$fitted.values <- NULL
    object$residuals <- NULL
    object$msg <- NULL
    ## allow a fixed routing model to be fit once
    object <- doParseRfit(object)
    object <- doRoutingFit(object, inverseFitOnly = TRUE)
    ## update parameters.
    ## the arguments in `...` may be intended for sma and/or routing
    dots <- list(...)
    ## first, store the SMA parameters to check whether they are changed.
    oldSmaPar <- coef(object, which = "sma", warn = FALSE)
    ## note that any NULL values in dots will delete parlist items
    ## (those parameters will then fall back to sim function defaults)
    object$parlist <- modifyList(object$parlist, dots)
    ## check whether the parameters are fully specified
    if (!isFullySpecified(object, which = "sma")) {
        object$U <- NULL
        return(object)
    }
    ## otherwise... SMA parameters are fully specified.
    ## need to run the SMA if parameters changed.
    newSmaPar <- coef(object, which = "sma")
    if (!identical(oldSmaPar, newSmaPar)) {
        RUNSMA <- TRUE
    }
    ## run SMA to generate U
    if (RUNSMA) {
        object$U <- predict(object, which = "sma")
    }
    ## handle routing
    routing <- object$routing
    if (!is.null(routing)) {
        ## take starting value of filter from data
        initX <- object$initX
        Xs_0 <- if (is.logical(initX)) 0 else initX
        if (is.null(object$initX)) {
            Q <- object$data[,"Q"]
            if (any(is.finite(Q[1:10]))) {
                object$initX <- min(Q[1:10], na.rm=TRUE)
            }
        }
        ## TODO: use initX!

        rfit <- object$rfit
        if (!is.null(rfit)) {
            ## fit the routing model.
            object <- doRoutingFit(object)
            ## it may have failed:
            if (!isFullySpecified(object))
                return(object)

        } else {
            ## just run routing simulation
            if (!isFullySpecified(object))
                return(object)

            ## TODO: take starting value of filter from data

            object$fitted.values <- predict(object, which = "routing")
#            object$residuals <-
#                (object$data[,"Q"] - object$fitted.values)
        }
    }
    ## absorbScale:
    ## move "steady state gain" (volume factor) from unit hydrograph
    ##      -- or model bias if no suitable UH -- into SMA
    ## (if possible; dispatched on SMA class)
    if (!and.rescale)
        return(object)
    alreadyRescaled <- FALSE
    if (!is.null(routing)) {
        rcoef <- coef(object, which = "routing")
        ssgFn <- paste("ssg", routing, sep = ".")
        if (exists(ssgFn, mode = "function")) {
            gain <- do.call(ssgFn, list(rcoef))
            if (!is.na(gain) && (abs(gain - 1) > 1e-6)) {
                ## UH needs to be normalised (unit volume)
                sobject <- absorbScale(object, gain = gain)
                if (!is.null(sobject)) {
                    object <- sobject
                }
                ## Even if we can't absorb scale into SMA:
                ## normalise UH parameters
                normFn <- paste("normalise", routing, sep = ".")
                if (exists(normFn, mode = "function")) {
                    normPars <- do.call(normFn, list(rcoef))
                    object$parlist <- modifyList(object$parlist,
                                                 as.list(normPars))
                } else {
                    stop("could not find function ", normFn,
                         " (works with ", ssgFn, ")")
                }
                ## note, we might not need to run the model again;
                ## assuming that absorbScale() on SMA has exactly compensated for
                ## normalise() on the routing parameters.
                object <- update(object, and.rescale = FALSE)
                ## NOTE - assuming that after normalise(), ssg.*() will return 1
                ##      - otherwise infinite recursion will result.
                alreadyRescaled <- TRUE
                ## TODO
            }
        }
    }

    ## TODO: i think should leave this out... ?

    if (!alreadyRescaled) {
        ## estimate scale by mass balance
        Q <- observed(object)
        X <- fitted(object)
        ok <- complete.cases(Q, X)
        gain <- (sum(Q[ok]) / sum(X[ok]))
        if (!is.finite(gain) || (gain == 0))
            return(object)
        sobject <- absorbScale(object, gain = gain)
        if (!is.null(sobject)) {
            object <- sobject
            ## and re-do the routing
            object <- update(object)
        }
    }

    return(object)
}

doParseRfit <-
    function(object)
{
    ## put 'rfit' into standard list form; AND
    ## calculate prefilter and add to 'rfit' specification
    routing <- object$routing
    rfit <- object$rfit
    if (is.null(routing))
        return(object)
    if (is.null(rfit))
        return(object)
    ## rfit was given
    if (identical(rfit, TRUE))
        rfit <- list()
    if (is.character(rfit))
        rfit <- list(rfit)
    if (!is.list(rfit))
        stop("unrecognised value of 'rfit'")

    ## TODO: set normalise = FALSE for armax/expuh if an absorbScale is available
    ## ?
    
    object$rfit <- rfit
    return(object)
}
