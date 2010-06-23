## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## Entries in hydromad.stats should all be functions with arguments (Q, X, ...)
## (and possible more: 'U', 'DATA', 'model').
## Expressions in the body of the function are allowed to be wrapped in .()
## in which case the fitting functions will cache the expression. It must
## only refer to 'Q' and/or 'DATA'.
## Note that hydromad.stats() inserts a function .() in these
## functions' environments which does nothing, therefore allowing them to be
## called directly even when the .() chunks have not been evaluated / cached.

## This is called once by fitting functions etc.
buildCachedObjectiveFun <-
    function(objective, model,
             DATA = observed(model, item = TRUE), Q = DATA[,"Q"])
{
    ## both 'Q' and 'DATA' can be referred to in .() expressions
    ## evaluate and replace .() expressions in body of objective
    if (inherits(objective, "formula")) {
        objective[[2]] <-
            eval(substitute(bquote(x), list(x = objective[[2]])))
    } else if (is.function(objective)) {
        body(objective) <- 
            eval(substitute(bquote(x), list(x = body(objective))))
    } else {
        stop("'objective' should be a function or formula, not a ",
             toString(class(objective)))
    }
    objective
}

## TODO: associate 'objective' with model rather than optimisation.
## then can build objective inside update.hydromad()...
## still - need to recalc stats for summary()?
## or could we store cached items in the model itself?
## --- store scalar model (fitted()) inside every hydromad object?

.defaultHydromadStats <- function()
    list("RMSE" = function(Q, X, ...) sqrt(mean((X - Q)^2, na.rm = TRUE)),
         "abs.err" = function(Q, X, ...) mean(abs(X - Q), na.rm = TRUE),
         "bias" = function(Q, X, ...) mean(X - Q, na.rm = TRUE),
         "rel.bias" = function(Q, X, ...) {
             ok <- complete.cases(X, Q)
             mean((X-Q)[ok]) / mean(Q[ok])
         },
         "r.squared" = function(Q, X, ...) {
             1 - fitStat(coredata(Q), coredata(X), ...)
         },
         "r.sq.boxcox" = function(Q, X, ...) {
             1 - .(buildObjectiveFun(Q, boxcox = TRUE))(Q, X, ...)
         },
         "r.sq.boxcox.old" = function(Q, X, ...) {
             1 - fitStat(Q, X, ..., trans = function(x) {
                 box.cox(x,
                         p = .(box.cox.powers(coredata(na.omit(Q)) +
                                              quantile(coredata(Q[Q>0]), 0.1, na.rm = TRUE)
                                              )$lambda),
                         start = .(quantile(coredata(Q[Q>0]), 0.1, na.rm = TRUE)))
             })
         },
         "r.sq.sqrt" = function(Q, X, ...) {
             1 - fitStat(Q, X, ..., trans = sqrt)
         },
         "r.sq.log" = function(Q, X, ...) {
             1 - fitStat(Q, X, ..., trans = function(x)
                         log(x + .(quantile(coredata(Q[which(Q>0)]), 0.1, na.rm = TRUE, names = FALSE))))
         },
         "r.sq.rank" = function(Q, X, ...) {
             1 - fitStat(Q, X, ..., trans = function(x) {
                 rank(round(log10(zapsmall(x, digits = 3)), digits = 2),
                      na.last = "keep")
             })
         },
         "r.sq.diff" = function(Q, X, ...) {
             1 - fitStat(diff(Q), diff(X), ...)
         },
         "r.sq.monthly" = function(Q, X, ...) {
             1 - .(buildObjectiveFun(Q, aggr.by = "months"))(Q, X, ...)
         },
         "r.sq.smooth7" = function(Q, X, ...) {
             1 - fitStat(Q, X, ..., trans = function(x)
                         simpleSmoothTs(x, width = 7, c = 2))
         },
         "r.sq.seasonal" = function(Q, X, ...) {
             1 - fitStat(Q, X, ref = .(ave(Q, months(time(Q)))), ...)
         },
         ## TODO: this is probably a bad idea
         "r.sq.vs.scalar" = function(Q, X, ..., model) {
             ref <- fitted(update(model, sma = "scalar"))
             1 - fitStat(Q, X, ref = ref, ...)
         },
         "r.sq.vs.tf" = function(Q, X, ..., DATA) {
             ref <-
                 .(fitted(hydromad(DATA, sma = "scalar", routing = "armax",
                                   rfit = list("sriv", order = c(2,1)))))
             1 - fitStat(Q, X, ref = ref, ...)
         },
         "r.sq.vs.tf.bc" = function(Q, X, ..., DATA) {
             objfun <- .({
                 ref <-
                     fitted(hydromad(DATA, sma = "scalar", routing = "armax",
                                     rfit = list("sriv", order = c(2,1))))
                 buildObjectiveFun(Q, ref = ref, boxcox = TRUE)
             })
             1 - objfun(Q, X, ...)
         },
#         "r.sq.vs.P" = function(Q, X, ..., P) {
#             rx <- filter(P, ar(Q, demean=FALSE)$ar, "r")
#             rx <- rx * (mean(Q, na.rm = TRUE) /
#                         mean(rx, na.rm = TRUE))
#             1 - fitStat(Q, X, ref = rx, ...)
#         },
         "persistence" = function(Q, X, ...) {
             1 - fitStat(Q, X, ref = lag(Q, -1), ...)
         },
         "events.medsums" = function(Q, X, ...) {
             objfun <-
                 .(buildEventObjectiveFun(Q, thresh = median(coredata(Q), na.rm = TRUE),
                                          mingap = 5, mindur = 5, all = TRUE, FUN = sum))
             objfun(Q, X, ...)
         },
         "events.medsums.vs.tf.bc" = function(Q, X, ..., DATA) {
             objfun <- .({
                 ref <-
                     fitted(hydromad(DATA, sma = "scalar", routing = "armax",
                                     rfit = list("sriv", order = c(2,1))))
                 buildEventObjectiveFun(Q, thresh = median(coredata(Q), na.rm = TRUE),
                                        mingap = 5, mindur = 5, all = TRUE, FUN = sum,
                                        boxcox = TRUE)
             })
             objfun(Q, X, ...)
         },
         "events.90sums" = function(Q, X, ...) {
             objfun <-
                 .(buildEventObjectiveFun(Q, thresh = quantile(coredata(Q), 0.9, na.rm = TRUE),
                                          mingap = 5, all = TRUE, FUN = sum))
             objfun(Q, X, ...)
         },
         "events.90sums.bc" = function(Q, X, ...) {
             objfun <-
                 .(buildEventObjectiveFun(Q, thresh = quantile(coredata(Q), 0.9, na.rm = TRUE),
                                          mingap = 5, all = TRUE, FUN = sum, boxcox = TRUE))
             objfun(Q, X, ...)
         },
         "events.90max" = function(Q, X, ...) {
             objfun <-
                 .(buildEventObjectiveFun(Q, thresh = quantile(coredata(Q), 0.9, na.rm = TRUE),
                                          mingap = 5, all = TRUE, FUN = max))
             objfun(Q, X, ...)
         },
         "events.90min" = function(Q, X, ...) {
             objfun <-
                 .(buildEventObjectiveFun(Q, thresh = quantile(coredata(Q), 0.9, na.rm = TRUE),
                                          mingap = 5, all = TRUE, FUN = min))
             objfun(Q, X, ...)
         },
         "ar1" = function(Q, X, ...) {
             cor(head(Q-X, -1), tail(Q-X, -1), use = "complete")
         },
         "X0" = function(Q, X, ...) cor(Q-X, X, use = "complete"),
         "X1" = function(Q, X, ...) cor(head(Q-X, -1), tail(X, -1), use = "complete"),
         "U1" = function(Q, X, ..., U) cor(head(Q-X, -1), tail(U, -1), use = "complete")
         )

## code below copied from lattice

hmadstat <- function(name, DATA = NULL)#, negate = FALSE)
{
    STATFUN <- .HydromadEnv$stats[[name]]
    if (!is.null(DATA))
        STATFUN <- buildCachedObjectiveFun(STATFUN, DATA = DATA)
#    if (negate)
#        return(function(...) - STATFUN(...))
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

    newfuns <- new[nm]
    
    ## ensure that everything being assigned is of the required type:
    ## functions of (Q, X, ...)
    if (!all(sapply(newfuns, is.function)))
        stop("items stored in hydromad.stats must be functions")
    if (!all(sapply(newfuns, function(x) "..." %in% names(formals(x)))))
        stop("items stored in hydromad.stats must be functions accepting '...'")
    ## insert a special function into its environment to allow direct evaluation:
    for (x in newfuns)
        assign(".", base::force, environment(x))

    ## this used to be

    ## modified <- updateList(retVal[nm], new[nm])
    ## .LatticeEnv$lattice.options[names(modified)] <- modified

    ## but then calling lattice.options(foo = NULL) had no effect
    ## because foo would be missing from modified.  So, we now do:

    updateList <- function (x, val) {
        if (is.null(x)) x <- list()
        utils::modifyList(x, val)
    }
    .HydromadEnv$stats <- updateList(old, newfuns)

    ## return changed entries invisibly
    invisible(retVal)
}
