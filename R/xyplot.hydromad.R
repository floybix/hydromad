## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


plot.hydromad <-
    function(x, y, ...)
{
    stop("There is no 'plot' method for 'hydromad' objects.",
         "Try 'xyplot', or 'plot(fitted(...))'")
}

xyplot.hydromad <-
    function(x, data = NULL,
             ...,
             all = FALSE,
             with.P = FALSE,
             screens = list(P = "rainfall", "streamflow"),
             type = list(P = "h", "l"),
             superpose = TRUE)
{
    stopifnot(is.null(data))
    tsdat <- cbind(modelled = fitted(x, all = all),
                   observed = observed(x, all = all))
    if (with.P)
        tsdat <- cbind(tsdat, P = observed(x, item = "P", all = all))
    ## NOTE superpose = FALSE is not really supported, because screens overrides it.
    foo <- xyplot(tsdat, ..., superpose = superpose,
                  screens = screens, type = type)
    foo$call <- sys.call(sys.parent())
    foo
}

xyplot.hydromad.runlist <-
    function(x, data = NULL,
             ...,
             all = FALSE,
             superpose = FALSE,
             with.P = FALSE,
             auto.key = TRUE,
             screens = list(P = "rainfall", "streamflow"),
             type = list(P = "h", "l"),
             layout = c(1, NA))
{
    stopifnot(is.null(data))
    if (superpose) {
        ## fitted models superposed, but rainfall still juxtaposed.
        ## include observed series from item 1 (assuming all are the same!)
        tsdat <- cbind(fitted(x, all = all),
                       observed = observed(x[[1]], all = all))
        if (with.P)
            tsdat <- cbind(tsdat, P = observed(x[[1]], item = "P", all = all))
        foo <- xyplot(tsdat, ..., superpose = superpose,
                      screens = screens, type = type, layout = layout)
    } else {
        ## fitted models juxtaposed, each with observed flow superposed
        if (!identical(auto.key, FALSE)) {
            if (isTRUE(auto.key)) auto.key <- list()
            if (is.null(auto.key$text))
                auto.key$text <- c("modelled", "observed")
        }
        foo <- xyplot.list(x, ..., all = all,
                           superpose = TRUE,
                           with.P = FALSE, type = type,
                           auto.key = auto.key,
                           layout = layout)
        if (with.P) {
            P <- observed(x[[1]], item = "P", all = all)
            rainPlot <-
                xyplot(P, type = "h", ...)
            foo <- c(foo, rainfall = rainPlot,
                     x.same = TRUE, y.same = NA, layout = layout)
        }
        
#        tsdat <- fitted(x, all = all)
#        foo <- xyplot(tsdat, ..., superpose = superpose,
#                      type = type,
#                      auto.key = list(text = c("modelled", "observed")),
#                      layout = layout)
#        ## include observed series from item 1 (assuming all are the same!)
#        foo <- foo + layer(panel.lines(obs),
#                           data = list(obs = observed(x[[1]], all = all)),
#                           style = 2, theme = trellis.par.get())
#        if (with.P) {
#            P <- observed(x[[1]], item = "P", all = all)
#            rainPlot <-
#                xyplot(P, type = "h",
#                       panel = function(...) {},
#                       ...) +
#                layer(panel.xyplot(...), style = 2, theme = trellis.par.get())
#            foo <- c(foo, rainfall = rainPlot,
#                     x.same = TRUE, y.same = NA, layout = layout)
#        }
    }
    foo$call <- sys.call(sys.parent())
    foo
}

qqmath.hydromad <-
    function(x, data = NULL,
             ...,
             all = FALSE,
             type = "l",
             auto.key = list(lines = TRUE, points = FALSE),
             f.value = ppoints(100), tails.n = 100)
{
    stopifnot(is.null(data))
    tsdat <- cbind(obs = observed(x, all = all),
                   mod = fitted(x, all = all))
    ## keep only common (corresponding) values
    tsdat[complete.cases(tsdat)==FALSE,] <- NA
    dat <- make.groups(observed = tsdat[,"obs"], modelled = tsdat[,"mod"])
    foo <- qqmath(~ data, groups = which, data = dat,
                  f.value = f.value, tails.n = tails.n,
                  auto.key = auto.key, type = type, ...)
    foo$call <- sys.call(sys.parent())
    foo
}

tsdiag.hydromad <- function(object, gof.lag, ...)
    stats:::tsdiag.Arima(object$uh, gof.lag = gof.lag, ...)

errormasscurve <- function(x, ...)
    UseMethod("errormasscurve")

errormasscurve.default <-
    function(x, ...)
{
    if (!is.atomic(x)) {
        x <- residuals(x)
        if (length(x) == 0) stop("could not get residuals() from 'x'")
    }
    bad <- is.na(x)
    x[bad] <- 0
    x[] <- cumsum(x)
    x[bad] <- NA
    foo <- xyplot(x, ...)
    foo$call <- sys.call(sys.parent())
    foo
}

