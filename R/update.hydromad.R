## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

update.hydromad <-
    function(object, ..., newdata = NULL,
             sma, routing, rfit, warmup, weights,
             and.rescale = TRUE)
{
    ## update timestamp
    object$last.updated <- Sys.time()

    if (!missing(sma))
        return(NextMethod("update"))
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
    if (!missing(routing))
        object$routing <- routing
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
        object$U <- predict(object, with_routing = FALSE)
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

            fnName <- paste(routing, "sim", sep = ".")
            force(get(fnName, mode = "function"))
            rcoef <- coef(object, which = "routing")
            rcall <- as.call(c(list(as.symbol(fnName),
                                    quote(object$U)),
                               as.list(rcoef)))
            object$fitted.values <- eval(rcall)
            object$residuals <-
                (object$data[,"Q"] - object$fitted.values)
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

absorbScale <- function(object, gain)
    UseMethod("absorbScale")

absorbScale.hydromad <- function(object, gain)
    return(NULL)

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

    ## TODO: set normalise = FALSE for expuh/uh if an absorbScale is available


    ## prefilter is relevant to these routing modules:
    if (routing %in% c("expuh", "uh")) {
        if (is.null(rfit$prefilter) &&
            isTRUE(hydromad.getOption("prefilter")))
        {
            rfit$prefilter <-
                makePrefilter(object$data, order = c(2,1))
        }
    }
    object$rfit <- rfit
    return(object)
}

doRoutingFit <-
    function(object, inverseFitOnly = FALSE)
{
    ## run 'rfit' with appropriate DATA
    routing <- object$routing
    rfit <- object$rfit
    if (is.null(routing))
        return(object)
    if (is.null(rfit))
        return(object)
    ## ok, 'rfit' given, fit the routing model.
    ## we should remove any existing routing parameters
    object$parlist <- as.list(coef(object, which = "sma", warn = FALSE))
    object$vcov.rfit <- NULL
    doRfit <- function(method = hydromad.getOption("rfit.method"), ...)
    {
        if (!is.character(method))
            stop("unrecognised 'method' in 'rfit' specification")
        isInverseMethod <- (any(grep("inverse", method)))
        if (isInverseMethod) {
            DATA <- object$data
        } else {
            if (inverseFitOnly)
                return(NULL)
            DATA <- cbind(U = object$U,
                          Q = object$data[,"Q"])
        }
        fnName <- paste(routing, method, "fit", sep = ".")
        force(get(fnName, mode = "function"))
        if (isInverseMethod && !isTRUE(hydromad.getOption("quiet")))
            message("== fitting routing by ", fnName, " method ==")
        fcall <- quote(RFIT(DATA, ...))
        fcall[[1]] <- as.symbol(fnName)
        ans <- eval(fcall)
        if (isInverseMethod && is.list(ans)) {
            ## these should be set to NULL anyway, but to be sure:
            ans$fitted.values <- NULL
            ans$residuals <- NULL
        }
        ans
    }
    rans <- do.call(doRfit, rfit)
    if (is.null(rans))
        return(object)
    if (inherits(rans, "try-error") ||
        inherits(rans, "error")) {
        object$msg <- rans
        return(object)
    }
    rcoef <- coef(rans)
    object$parlist <- modifyList(object$parlist, as.list(rcoef))
    object$fitted.values <- fitted(rans)
    object$residuals <- residuals(rans)
    tmp <- try(vcov(rans), silent = TRUE)
    if (inherits(tmp, "try-error"))
        tmp <- rans$vcov
    object$vcov.rfit <- tmp
    ## delete the fit specification, as it has done its work now
    object$used.rfit <- rfit
    object$rfit <- NULL
    object
}
