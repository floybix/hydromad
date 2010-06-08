## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

.defaultHydromadStats <- function()
    list("RMSE" = function(Q, X, ...) sqrt(mean((X - Q)^2, na.rm = TRUE)),
         "abs.err" = function(Q, X, ...) mean(abs(X - Q), na.rm = TRUE),
         "bias" = function(Q, X, ...) mean(X - Q, na.rm = TRUE),
         "rel.bias" = function(Q, X, ...) {
             ok <- complete.cases(X, Q)
             mean((X-Q)[ok]) / mean(Q[ok])
         },
         "r.squared" = function(Q, X, ...) fitStat(Q, X),
         "r.sq.sqrt" = function(Q, X, ...) fitStat(Q, X, trans = sqrt),
         "r.sq.log" = function(Q, X, ...) fitStat(Q, X, trans = log),
         "r.sq.rank" = function(Q, X, ...) {
             fitStat(Q, X, trans = function(x) {
                 rank(round(log10(zapsmall(x, digits = 3)), digits = 2),
                      na.last = "keep")
             })
         },
         "r.sq.diff" = function(Q, X, ...) {
             fitStat(diff(Q), diff(X))
         },
         "r.sq.monthly" = function(Q, X, ...) {
             tsFitStat(Q, X, aggr = list(by = cut(time(Q), "months"), FUN = sum))
         },
         "r.sq.smooth7" = function(Q, X, ...) {
             tsFitStat(Q, X, trans = function(x) simpleSmoothTs(x, width = 7, c = 2))
         },
         "r.sq.seasonal" = function(Q, X, ...) {
             tsFitStat(Q, X, ref = ave(Q, months(time(Q))))
         },
         "r.sq.vs.P" = function(Q, X, ..., P) {
             rx <- filter(P, ar(Q, demean=FALSE)$ar, "r")
             rx <- rx * (mean(Q, na.rm = TRUE) /
                         mean(rx, na.rm = TRUE))
             fitStat(Q, X, ref = rx)
         },
         "persistence" = function(Q, X, ...) {
             tsFitStat(Q, X, ref = lag(Q, -1))
         },
         "events.medsums" = function(Q, X, ...) {
             tsFitStat(Q, X, 
                       events = list(thresh = median(coredata(Q), na.rm = TRUE),
                       mingap = 5, mindur = 5, all = TRUE, FUN = sum))
         },
         "events.90sums" = function(Q, X, ...) {
             tsFitStat(Q, X, 
                       events = list(thresh = quantile(coredata(Q), 0.9, na.rm = TRUE),
                       mingap = 5, all = TRUE, FUN = sum))
         },
         "events.90max" = function(Q, X, ...) {
             tsFitStat(Q, X, 
                       events = list(thresh = quantile(coredata(Q), 0.9, na.rm = TRUE),
                       mingap = 5, FUN = max))
         },
         "events.90min" = function(Q, X, ...) {
             tsFitStat(Q, X, 
                       events = list(thresh = quantile(coredata(Q), 0.9, na.rm = TRUE),
                       below = TRUE, mindur = 5, FUN = min))
         },
         "ar1" = function(Q, X, ...) {
             cor(head(Q-X, -1), tail(Q-X, -1), use = "complete")
         },
         "X0" = function(Q, X, ...) cor(Q-X, X, use = "complete"),
         "X1" = function(Q, X, ...) cor(head(Q-X, -1), tail(X, -1), use = "complete"),
         "U1" = function(Q, X, ..., U) cor(head(Q-X, -1), tail(U, -1), use = "complete")
         )

## code below copied from lattice

hmadstat <- function(name, negate = FALSE)
{
    STATFUN <- .HydromadEnv$stats[[name]]
    if (negate)
        return(function(...) - STATFUN(...))
    STATFUN
}

hydromad.stats <- function(...)
{
    ## this would have been really simple if only form allowed were
    ## lattice.options("foo", "bar") and
    ## lattice.options(foo=1, bar=2). But it could also be
    ## lattice.options(foo=1, "bar"), which makes some juggling necessary

    new <- list(...)
    if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) new <- new[[1]]
    old <- .HydromadEnv$stats

    ## if no args supplied, returns full options list
    if (length(new) == 0) return(old)

    nm <- names(new)
    if (is.null(nm)) return(old[unlist(new)]) ## typically getting options, not setting
    isNamed <- nm != "" ## typically all named when setting, but could have mix
    if (any(!isNamed)) nm[!isNamed] <- unlist(new[!isNamed])

    ## so now everything has non-"" names, but only the isNamed ones should be set
    ## everything should be returned, however

    retVal <- old[nm]
    names(retVal) <- nm
    nm <- nm[isNamed]

    ## this used to be

    ## modified <- updateList(retVal[nm], new[nm])
    ## .LatticeEnv$lattice.options[names(modified)] <- modified

    ## but then calling lattice.options(foo = NULL) had no effect
    ## because foo would be missing from modified.  So, we now do:

    updateList <- function (x, val) {
        if (is.null(x)) x <- list()
        utils::modifyList(x, val)
    }
    .HydromadEnv$stats <- updateList(old, new[nm])

    ## return changed entries invisibly
    invisible(retVal)
}
