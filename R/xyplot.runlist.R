## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## This plots fitted vs observed only.
## To plot residuals, just call xyplot(residuals())
xyplot.runlist <-
    function(x, data = NULL,
             ...,
             all = FALSE,
             superpose = FALSE,
             auto.key = list(points = FALSE, lines = TRUE),
             x.same = TRUE, y.same = NA, layout = c(1, NA))
{
    if (superpose) {
        ## include observed series from item 1 (assuming all are the same!)
        tsdat <- cbind(fitted(x, all = all),
                       observed = observed(x[[1]], all = all))
        foo <- xyplot(tsdat, ..., superpose = TRUE,
                      auto.key = auto.key, layout = layout)
    } else {
        foo <- xyplot.list(x, ..., auto.key = auto.key,
                           x.same = x.same, y.same = y.same, layout = layout)
    }
    foo$call <- sys.call(sys.parent())
    foo
}

## Handles either fitted vs observed, or residuals.
qqmath.runlist <-
    function(x, data = NULL,
             ...,
             residuals = FALSE,
             all = FALSE,
             superpose = FALSE,
             f.value = ppoints(100), tails.n = 100,
             type = "l",
             auto.key = list(lines = TRUE, points = FALSE))
{
    if (superpose) {
        if (residuals) {
            tsdat <- residuals(x, all = all)
        } else {
            ## include observed series from item 1 (assuming all are the same!)
            tsdat <- cbind(fitted(x, all = all),
                           observed = observed(x[[1]], all = all))
        }
        dat <- do.call("make.groups", as.data.frame(tsdat))
        foo <- qqmath(~ data, groups = which, data = dat,
                      f.value = f.value, tails.n = tails.n,
                      type = type, auto.key = auto.key)
    } else {
        if (residuals) {
            FUN <- function(x, ...) qqmath(~ residuals(x), ...)
        } else {
            FUN <- qqmath
        }
        foo <- xyplot.list(x, ..., FUN = FUN,
                           f.value = f.value, tails.n = tails.n,
                           type = type, auto.key = auto.key)
    }
    foo$call <- sys.call(sys.parent())
    foo
}
