## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## aspects:
## * aggregation (regular / events)
## * reference model
## * transformation

buildObjectiveFun <-
    function(Q, groups = NULL,
             #aggr.by = NULL,
             ..., FUN = sum,
             ref = NULL,
             boxcox = FALSE, start = NULL)
{
    if (is.null(aggr.by)) {
        doaggr <- identity
    } else {
        doaggr <- function(x)
            eventapply(x, groups, FUN = FUN)
    }
    aggrQ <- doaggr(Q)
    aggrRef <- NULL
    if (!is.null(ref))
        aggrRef <- doaggr(ref)
    if (boxcox) {
        if (is.null(start))
            start <-
                quantile(coredata(aggrQ[which(aggrQ > 0)]), 0.1,
                         na.rm = TRUE, names = FALSE)
        lambda <- box.cox.powers(coredata(na.omit(aggrQ)) + start)$lambda
        function(Q, X, ...) {
            fitStat(aggrQ, doaggr(X), ref = aggrRef, ...,
                    trans = function(x) box.cox(x, lambda, start = start))
        }
    } else {
        function(Q, X, ...) {
            fitStat(aggrQ, doaggr(X), ref = aggrRef, ...)
        }
    }
}
