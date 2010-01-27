## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


xyplot.ihacres <-
xyplot.tf <-
    function(x, data,
             coerce = byDays, trans = NULL,
             superpose = TRUE,
             ...)
{
    if (!missing(data) && !is.null(data))
        warning("'data' argument ignored.")
    tsdat <- cbind(obs=observed(x), mod=fitted(x))
    tsdat <- coerce(tsdat)
    ## TODO: trans?
    foo <- xyplot(tsdat, superpose = superpose, ...)
    foo$call <- match.call()
    foo
}

errormasscurve <- function(x, ...)
    UseMethod("errormasscurve")

errormasscurve.default <-
    function(x, coerce=as.ts, ...)
{
    if (!is.atomic(x)) {
        x <- residuals(x)
        if (length(x) == 0) stop("could not get residuals() from 'x'")
    }
    x <- coerce(x)
    x[is.na(x)] <- 0
    foo <- xyplot(cumsum(x), ...)
    foo$call <- match.call() #sys.call(sys.parent())
    foo
}

qqmath.ihacres <-
    function(x, data,
             coerce = byDays, trans = NULL,
             type = c("p","l"),
             auto.key = TRUE,
             subset = complete.cases(obsmod), ...)
{
    if (!missing(data) && !is.null(data))
        warning("'data' argument ignored.")
    obsmod <- cbind(obs=observed(x), mod=fitted(x))
    obsmod <- coerce(obsmod)
    transfn <- eval(trans)
    if (is.character(transfn)) transfn <- get(trans)
    if (!is.null(trans))
        obsmod <- transfn(obsmod)
    dat <- make.groups(observed=obsmod[,"obs"], modelled=obsmod[,"mod"])
    foo <- qqmath(~ data, groups=which, data=dat, auto.key=auto.key,
                  type=type, subset=subset, ...)
    foo$call <- sys.call(sys.parent())
    foo
}

tsdiag.ihacres <- function(object, gof.lag, ...)
    stats:::tsdiag.Arima(object$uh, gof.lag=gof.lag, ...)

