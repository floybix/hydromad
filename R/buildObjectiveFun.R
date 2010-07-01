## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## aspects:
## * aggregation groups (regular / events) & aggregation function
## * reference model
## * data transformation

buildObjectiveFun <-
    function(Q, groups = NULL, FUN = sum, ...,
             ref = NULL, boxcox = FALSE, start = NULL)
{
    doaggr <- identity
    if (!is.null(groups)) {
        argsForFUN <- list(...)
        fullFUN <- function(...)
            do.call(FUN, modifyList(list(...), argsForFUN))
        doaggr <- function(x)
            eventapply(x, groups, FUN = fullFUN)
    }
    aggrQ <- doaggr(Q)
    aggrRef <- NULL
    if (!is.null(ref))
        aggrRef <- doaggr(ref)
    if (!identical(boxcox, FALSE) && (length(aggrQ) > 1)) {
        coreaggrQ <- coredata(na.omit(aggrQ))
        if (is.null(start))
            start <-
                quantile(coreaggrQ[coreaggrQ > 0], 0.1, names = FALSE)
        if (isTRUE(boxcox)) {
            lambda <- box.cox.powers(coreaggrQ + start)$lambda
        } else {
            stopifnot(is.numeric(boxcox))
            lambda <- boxcox
        }
        function(Q, X, ...) {
            nseStat(aggrQ, doaggr(X), ref = aggrRef, ...,
                    trans = function(x) box.cox(x, lambda, start = start))
        }
    } else {
        function(Q, X, ...) {
            nseStat(aggrQ, doaggr(X), ref = aggrRef, ...)
        }
    }
}
