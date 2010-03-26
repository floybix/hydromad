## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


xyplot.hydromad <-
    function(x, data = NULL,
             coerce = byDays, all = FALSE,
             superpose = TRUE,
             with.P = FALSE,
             screens = list(P = "rainfall", "streamflow"),
             ...)
{
    stopifnot(is.null(data))
    tsdat <- cbind(obs = observed(x, all = all), mod = fitted(x, all = all))
    if (with.P)
        tsdat <- cbind(tsdat, P = x$data[,"P"])
    tsdat <- coerce(tsdat)
    foo <- xyplot(tsdat, superpose = superpose, screens = screens, ...)
    foo$call <- sys.call(sys.parent())
    foo
}

errormasscurve <- function(x, ...)
    UseMethod("errormasscurve")

errormasscurve.default <-
    function(x, coerce = as.ts, ...)
{
    if (!is.atomic(x)) {
        x <- residuals(x)
        if (length(x) == 0) stop("could not get residuals() from 'x'")
    }
    x <- coerce(x)
    x[is.na(x)] <- 0
    x[] <- cumsum(x)
    foo <- xyplot(x, ...)
    foo$call <- sys.call(sys.parent())
    foo
}

qqmath.hydromad <-
    function(x, data,
             coerce = byDays, trans = NULL,
             type = c("p","l"),
             auto.key = list(lines = TRUE, points = FALSE),
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

tsdiag.hydromad <- function(object, gof.lag, ...)
    stats:::tsdiag.Arima(object$uh, gof.lag=gof.lag, ...)

