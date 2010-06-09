

optimtrace <- function(object, ...)
    UseMethod("optimtrace")

optimtrace.default <- function(object, ..., model = object)
{
    stopifnot(is.list(object))
    stopifnot(length(object$fit.result) > 0)

    getObjSeq(object$fit.result, ..., model = model)
}

getObjSeq <- function(object, ...)
    UseMethod("getObjSeq")

## method for optim() results produced by hydromad
getObjSeq.default <- function(object, ..., raw = FALSE)
{
    ppp <- object$objseq
    ans <- zoo(ppp, 1:NROW(ppp))
    if (raw == FALSE) ans <- cummin(na.locf(ans))
    ans
}

## TODO: could re-calculate objective using stored population
getObjSeq.SCEoptim <- function(object, ..., raw = FALSE)
{
    funevals <- object$counts
    ppp <- t(object$POP.FIT.ALL)
    if (raw == FALSE) {
        ## best of the chains at each step
        ppp <- do.call("pmin", as.data.frame(ppp))
    }
    #ts(ppp, end = funevals, deltat = funevals / nrow(ppp))
    zoo(ppp, seq(0, funevals, length = NROW(ppp)+1)[-1])
}

## TODO: could re-calculate objective using stored population
getObjSeq.DEoptim <- function(object, ...)
{
    funevals <- object$optim$nfeval
    ppp <- object$member$bestvalit
    zoo(ppp, seq(0, funevals, length = NROW(ppp)+1)[-1])
}

getObjSeq.dream <-
    function(object, ..., raw = FALSE, objective = NULL, model = NULL)
{
    funevals <- object$fun.evals
    ppp <- - object$hist.logp ## negative log likelihood
    if (is.null(objective)) {
        if (raw == FALSE) {
            ## best of the chains at each step
            ppp <- do.call("pmin", as.data.frame(ppp))
            ## best result so far at each step
            ppp <- cummin(na.locf(ppp))
        }
    } else {
        ## calculate corresponding objective function values over time.
        stopifnot(!is.null(model))
        ## could do it in a straightforward way by evaluating model at
        ## each time step and for each chain, but would take a long time!
        ## so instead, only re-calculate objective when likelihood improves.
        bestLik <- do.call("pmin", as.data.frame(ppp))
        improved <- c(TRUE, diff(cummin(na.locf(bestLik))) != 0)
        ## calculate objective for time steps with improved likelihood
        improv.fval <- sapply(which(improved), function(i) {
            ichain <- which.min(ppp[i,])
            ipars <- object$Sequences[[ichain]][i,]
            objFunVal(update(model, newpars = ipars),
                      objective = objective)
        })
        fval <- rep(NA_real_, NROW(ppp))
        fval[improved] <- improv.fval
        fval <- na.locf(fval)
        ppp <- fval
    }
    zoo(ppp, seq(0, funevals, length = NROW(ppp)+1)[-1])
}
