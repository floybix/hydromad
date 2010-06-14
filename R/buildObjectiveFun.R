## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## aspects:
## * aggregation (regular / events)
## * reference model
## * transformation

buildObjectiveFun <-
    function(Q, aggr.by = NULL, ..., FUN = sum,
             ref = NULL,
             boxcox = FALSE, start = NULL)
{
    if (is.null(aggr.by)) {
        doaggr <- identity
    } else if (identical(aggr.by, "events")) {
        events <- eventseq(Q, ...)
        doaggr <- function(x) eventapply(x, events, FUN = FUN)
    } else {
        groups <- cut(time(Q), aggr.by)
        ## convert to character
        ## TODO - this should not be necessary - to be fixed in zoo?
        #groups <- levels(groups)[groups]
        groups <- unclass(groups)
#        if (inherits(time(Q), "Date"))
#            groups <- as.Date(levels(groups))[groups]
#        else if (inherits(time(Q), "POSIXct"))
#            groups <- as.POSIXct(levels(groups))[groups]
#        else
#            groups <- (time(Q)[!duplicated(groups)])[groups]
        doaggr <- function(x) aggregate(x, groups, FUN = FUN)
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

buildEventObjectiveFun <-
    function(Q, thresh = 0, ..., FUN = sum,
             ref = NULL,
             boxcox = FALSE, start = NULL)
{
    buildObjectiveFun(Q, aggr.by = "events",
                      thresh = thresh, ..., FUN = FUN,
                      ref = ref, boxcox = boxcox, start = start)
}
