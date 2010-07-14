

event.xyplot <- function(x, ...)
    UseMethod("event.xyplot")

event.xyplot.formula <-
    function(x, data = NULL, events, ...,
             response = NULL, eFUN = sum,
             as.table = TRUE, layout = NULL,
             auto.key = NCOL(response) > 1,
             abline = list(h = 0),
             xlab = "event-aggregated covariate values",
             ylab = "event-aggregated response")
{
    x <- as.formula(x)
    ## define a special function to aggregate events, for use in formula
    e <- function(z, FUN = sum, ...) eventapply(z, events, FUN = FUN, ...)
    ## make it known to the formula:
    environment(x) <- environment()
    ## another convenience function for formula:
    first <- function(x) head(x, 1)
    ## extract covariates from formula
    tt <- terms(x, data = data, keep.order = TRUE)
    mf <- model.frame(tt, data)
    covars <- model.matrix(tt, mf)
    ## remove automatic intercept
    covars <- covars[, -1, drop = FA:SE]
    if (is.null(response)) {
        response <- model.response(mf)
    } else {
        response <- e(response, eFUN)
    }
    ystack <- stack(as.data.frame(response))
    
    foo <-
        xyplot.list(as.data.frame(covars), data = ystack,
                    FUN = function(VAR, ...)
                    xyplot(values ~ rep(VAR, NCOL(response)), groups = ind, ...),
                    default.scales = list(x = list(relation = "free")),
                    ...,
                    as.table = as.table, layout = layout,
                    auto.key = auto.key, 
                    abline = abline,
                    xlab = xlab, ylab = ylab)
    foo$call <- sys.call(sys.parent())
    foo
}

event.xyplot.hydromad.runlist <-
event.xyplot.hydromad <-
    function(x, events,
             formula = ~ e(time(data),first) + log10(e(Q,sum)+1) +
                         sqrt(e(P,max)) + log10(e(lag(Q,-2),first)+1),
             extract = function(x) residuals(x, boxcox = FALSE),
             ...,
             panel = function(x ,y, ..., se) {
                 panel.rug(x = x, alpha = 0.5)
                 panel.abline(h = 0)
                 library(mgcv)
                 panel.smoother(x, y, y ~ s(x), method = "gam", se = se)
                 panel.xyplot(x, y, ...)
             },
             pch = ".", se = FALSE,
             ylab = "residual flow sums in event windows (mm)",
             data = NULL)
{
    if (!is.null(data)) warning("'data' ignored")
    if (inherits(x, "runlist")) {
        data <- observed(x[[1]], select = TRUE)
    } else {
        data <- observed(x, select = TRUE)
        data$U <- fitted(x, U = TRUE)
    }
    data$julian <- julian(time(data))
    response <- extract(x)
    foo <- 
        event.xyplot(formula, data = data,
                     events = events,
                     response = response,
                     ...,
                     panel = panel,
                     pch = pch, se = se,
                     ylab = ylab)
    foo$call <- sys.call(sys.parent())
    foo
}
