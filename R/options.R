## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

.defaultHydromadOptions <- function()
    list(
         sma = NULL,
         routing = "uh",
         cwi = cwi.ranges(),
         cmd = cmd.ranges(),
         sacramento = sacramento.ranges(),
         bucket = bucket.ranges(),
         rfit.method = "sriv",
         inverse.fit.method = "sriv",
         order = c(n = 1, m = 0),
         delay = NA,
         max.delay = 10,
         warmup = 100,
         normalise = TRUE,
         fit.samples = 64,
         objective = ~ -2 * fitStat(Q, X, trans = sqrt) + abs(fitBias(Q, X)),
         prefilter = TRUE,
         riv.noise.order = NULL,
         sriv.iterations = 16,
         sriv.epsilon = 1e-3,
         inverse.iterations = 30,
         inverse.tolerance = 1e-5,
         inverse.epsilon = 0.001,
         optim.method = "BFGS",
         optim.control = list(fnscale = 1, reltol = 1e-5, maxit = 150,
                trace = if (interactive()) 4 else 0, REPORT = 4),
         uhParTrans =
             alist(tau = list(trans = log, inverse = exp),
                   v = list(trans = qlogis, inverse = plogis),
                   #v = list(trans = function(x) qlogis(max(x, 0.001)),
                   #         inverse = plogis),
                   #v = list(trans = log, inverse = exp),
                   lambda = list(trans = function(x) qlogis(-x),
                                inverse = function(x) -plogis(x)),
                   loss = list(inverse = function(x) max(x, 0))),
         sim.epsilon = 1e-5,
         quiet = FALSE,
         trace = FALSE,
         catch.errors = TRUE,
         catch.errors.optim = TRUE,
         pure.R.code = FALSE
         )

## code below copied from lattice

hydromad.getOption <- function(name)
{
    .HydromadEnv$options[[name]]
}

hydromad.options <- function(...)
{
    ## this would have been really simple if only form allowed were
    ## lattice.options("foo", "bar") and
    ## lattice.options(foo=1, bar=2). But it could also be
    ## lattice.options(foo=1, "bar"), which makes some juggling necessary

    new <- list(...)
    if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) new <- new[[1]]
    old <- .HydromadEnv$options

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
    .HydromadEnv$options <- updateList(old, new[nm])

    ## return changed entries invisibly
    invisible(retVal)
}
