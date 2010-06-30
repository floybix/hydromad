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




eventseq_q90 <-
    function(Q, continue = !all, all = FALSE)
{
    q90 <- quantile(coredata(Q), 0.9, na.rm = TRUE)
    ev <- eventseq(Q, thresh = q90, indur = 4,
                   mindur = 5, mingap = 5, continue = continue, all = all)
    ev
}



## This is called once by fitting functions etc.
buildCachedObjectiveFun <-
    function(objective, model,
             DATA = observed(model, item = TRUE), Q = DATA[,"Q"])
{
    ## both 'Q' and 'DATA' can be referred to in .() expressions
    ## evaluate and replace .() expressions in body of objective
    if (inherits(objective, "formula")) {
        if (length(objective) > 2)
            warning("left hand side of formula ignored")
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

.defaultHydromadStats <- function()
{
    ## keep R CMD check happy:
    . <- function(x) x
    ## default set of stats
    list("bias" = function(Q, X, ...) mean(X - Q, na.rm = TRUE),
         "rel.bias" = function(Q, X, ...) {
             ok <- complete.cases(X, Q)
             mean((X-Q)[ok]) / mean(Q[ok])
         },
         "abs.err" = function(Q, X, ...) mean(abs(X - Q), na.rm = TRUE),
         "RMSE" = function(Q, X, ...) sqrt(mean((X - Q)^2, na.rm = TRUE)),
         "r.squared" = function(Q, X, ...) {
             1 - fitStat(coredata(Q), coredata(X), ...)
         },
         "r.sq.sqrt" = function(Q, X, ...) {
             1 - fitStat(Q, X, ..., trans = sqrt)
         },
         "r.sq.log" = function(Q, X, ...) {
             1 - fitStat(Q, X, ..., trans = function(x)
                         log(x + .(quantile(coredata(subset(Q, Q>0)), 0.1, na.rm = TRUE, names = FALSE))))
         },
         "r.sq.boxcox" = function(Q, X, ...) {
             1 - .(buildObjectiveFun(Q, boxcox = TRUE))(Q, X, ...)
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
             1 - .(buildObjectiveFun(Q, groups = cut(time(Q), "months")))(Q, X, ...)
         },
         "r.sq.smooth5" = function(Q, X, ...) {
             1 - fitStat(Q, X, ..., trans = function(x)
                         simpleSmoothTs(x, width = 5, c = 2))
         },
         "r.sq.seasonal" = function(Q, X, ...) {
             1 - fitStat(Q, X,
                         ref = .(ave(Q, months(time(Q)), FUN = function(x) mean(x, na.rm = TRUE))),
                         ...)
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
         "persistence" = function(Q, X, ...) {
             1 - fitStat(Q, X, ref = lag(Q, -1), ...)
         },
         "persistence.bc" = function(Q, X, ...) {
             objfun <- .({
                 ref <- lag(Q, -1)
                 buildObjectiveFun(Q, ref = ref, boxcox = TRUE)
             })
             1 - objfun(Q, X, ...)
         },
         
         "e.rain5" = function(Q, X, ..., DATA) {
             objfun <- .({
                 ev <- eventseq(DATA$P, thresh = 5, inthresh = 1,
                                indur = 4, continue = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = sum)
             })
             objfun(Q, X, ...)
         },
         "e.rain5.log" = function(Q, X, ..., DATA) {
             objfun <- .({
                 ev <- eventseq(DATA$P, thresh = 5, inthresh = 1,
                                indur = 4, continue = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = sum, boxcox = 0)
             })
             objfun(Q, X, ...)
         },
         "e.rain5.bc" = function(Q, X, ..., DATA) {
             objfun <- .({
                 ev <- eventseq(DATA$P, thresh = 5, inthresh = 1,
                                indur = 4, continue = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = sum, boxcox = TRUE)
             })
             objfun(Q, X, ...)
         },
         
         #= flowq90_indur4_mindur5_mingap5

         "e.q90" = function(Q, X, ...) {
             objfun <-.({
                 q90 <- quantile(coredata(Q), 0.9, na.rm = TRUE)
                 ev <- eventseq(Q, thresh = q90, indur = 4,
                                mindur = 5, mingap = 5, continue = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = sum)
             })
             objfun(Q, X, ...)
         },
         "e.q90.log" = function(Q, X, ...) {
             objfun <-.({
                 q90 <- quantile(coredata(Q), 0.9, na.rm = TRUE)
                 ev <- eventseq(Q, thresh = q90, indur = 4,
                                mindur = 5, mingap = 5, continue = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = sum, boxcox = 0)
             })
             objfun(Q, X, ...)
         },
         "e.q90.bc" = function(Q, X, ...) {
             objfun <-.({
                 q90 <- quantile(coredata(Q), 0.9, na.rm = TRUE)
                 ev <- eventseq(Q, thresh = q90, indur = 4,
                                mindur = 5, mingap = 5, continue = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = sum, boxcox = TRUE)
             })
             objfun(Q, X, ...)
         },

         "e.q90.all" = function(Q, X, ...) {
             objfun <-.({
                 q90 <- quantile(coredata(Q), 0.9, na.rm = TRUE)
                 ev <- eventseq(Q, thresh = q90, indur = 4,
                                mindur = 5, mingap = 5, all = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = sum)
             })
             objfun(Q, X, ...)
         },
         "e.q90.all.log" = function(Q, X, ...) {
             objfun <-.({
                 q90 <- quantile(coredata(Q), 0.9, na.rm = TRUE)
                 ev <- eventseq(Q, thresh = q90, indur = 4,
                                mindur = 5, mingap = 5, all = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = sum, boxcox = 0)
             })
             objfun(Q, X, ...)
         },
         "e.q90.all.bc" = function(Q, X, ...) {
             objfun <-.({
                 q90 <- quantile(coredata(Q), 0.9, na.rm = TRUE)
                 ev <- eventseq(Q, thresh = q90, indur = 4,
                                mindur = 5, mingap = 5, all = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = sum, boxcox = TRUE)
             })
             objfun(Q, X, ...)
         },

         "e.q90.min" = function(Q, X, ...) {
             objfun <-.({
                 q90 <- quantile(coredata(Q), 0.9, na.rm = TRUE)
                 ev <- eventseq(Q, thresh = q90, indur = 4,
                                mindur = 5, mingap = 5, continue = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = min, na.rm = TRUE)
             })
             objfun(Q, X, ...)
         },
         "e.q90.min.log" = function(Q, X, ...) {
             objfun <-.({
                 q90 <- quantile(coredata(Q), 0.9, na.rm = TRUE)
                 ev <- eventseq(Q, thresh = q90, indur = 4,
                                mindur = 5, mingap = 5, continue = TRUE)
                 buildObjectiveFun(Q, groups = ev, FUN = min, na.rm = TRUE,
                                   boxcox = 0)
             })
             objfun(Q, X, ...)
         },
         
         "ar1" = function(Q, X, ...) {
             cor(head(Q-X, -1), tail(Q-X, -1), use = "complete")
         },
         "X0" = function(Q, X, ...) cor(Q-X, X, use = "complete"),
         "X1" = function(Q, X, ...) cor(head(Q-X, -1), tail(X, -1), use = "complete"),
         "U1" = function(Q, X, ..., U) cor(head(Q-X, -1), tail(U, -1), use = "complete")
         )
}

hmadstat <- function(name, DATA = NULL)
{
    STATFUN <- .HydromadEnv$stats[[name]]
    if (!is.null(DATA))
        STATFUN <- buildCachedObjectiveFun(STATFUN, DATA = DATA)
    STATFUN
}

## this function based on lattice::lattice.options, but with some extra bits

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
