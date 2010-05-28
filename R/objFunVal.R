## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


statDefns <- function()
{
    list(
         "bias" = ~ mean(X - Q, na.rm = TRUE),
         "rel.bias" =
         ~ { ok <- complete.cases(X,Q); mean((X-Q)[ok]) / mean(Q[ok]) },
         "r.squared" = ~ fitStat(Q, X),
         "r.sq.sqrt" = ~ fitStat(Q, X, trans = sqrt),
         "r.sq.log" = ~ fitStat(Q, X, trans = log, offset = TRUE),
         "r.sq.diff" = ~ fitStat(diff(Q), diff(X)),
         "r.sq.rank" =
         ~ fitStat(Q, X, trans = function(x)
                   rank(round(log10(zapsmall(x, digits = 3)), digits = 2),
                        na.last = "keep")),
         "r.sq.monthly" =
         ~ tsFitStat(Q, X, aggr = list(by = cut(time(Q), "months"), FUN = sum)),
         "r.sq.smooth7" =
         ~ tsFitStat(Q, X, trans = function(x) simpleSmoothTs(x, width = 7, c = 2)),
         "r.sq.seasonal" =
         ~ tsFitStat(Q, X, ref = ave(Q, months(time(Q)))),
         "r.sq.vs.P" =
         ~ fitStat(Q, X, ref = { rx <- filter(P, ar(Q, demean=FALSE)$ar, "r");
                                 rx * mean(Q, na.rm = TRUE) /
                                     mean(rx, na.rm = TRUE) }),
         "persistence" = ~ tsFitStat(Q, X, ref = lag(Q, -1)),
         "events.medsums" =
         ~ tsFitStat(Q, X, events = list(thresh = median(Q, na.rm = TRUE),
                           mingap = 5, mindur = 5, all = TRUE, FUN = sum)),
         "events.90sums" =
         ~ tsFitStat(Q, X, events = list(thresh = quantile(Q, 0.9, na.rm = TRUE),
                           mingap = 5, all = TRUE, FUN = sum)),
         "events.90max" =
         ~ tsFitStat(Q, X, events = list(thresh = quantile(Q, 0.9, na.rm = TRUE),
                           mingap = 5, FUN = max)),
         "events.90min" =
         ~ tsFitStat(Q, X, events = list(thresh = quantile(Q, 0.9, na.rm = TRUE),
                           below = TRUE, mindur = 5, FUN = min)),
         "abs.err" = ~ mean(abs(X - Q), na.rm = TRUE),
         "RMSE" = ~ sqrt(mean((Q - X)^2, na.rm = TRUE)),
         "ar1" = ~ cor(head(Q-X, -1), tail(Q-X, -1), use = "complete"),
         "X0" = ~ cor(Q-X, X, use = "complete"),
         "X1" = ~ cor(head(Q-X, -1), tail(X, -1), use = "complete"),
         "U1" = ~ cor(head(Q-X, -1), tail(U, -1), use = "complete")
         )
}

objFunVal <- function(x, objective, ...)
    UseMethod("objFunVal")

objFunVal.default <-
    function(x, objective = hydromad.getOption("objective"),
             ..., nan.ok = FALSE)
{
    stopifnot(is.numeric(x) || is.data.frame(x))
    stopifnot(length(colnames(x)) > 0)
    ## these can be referred to in `objective`
    X <- x[,"X"]
    delayedAssign("Q", x[,"Q"])
    delayedAssign("U", x[,"U"])
    delayedAssign("P", x[,"P"])
    objFunVal1 <- function(obj)
    {
        if (inherits(obj, "formula"))
            obj <- obj[[2]]
        stopifnot(is.language(obj))
        val <- eval(obj)
        if (is.nan(val)) {
            if (identical(nan.ok, "warn"))
                warning("objective function returned NaN")
            else if (isTRUE(nan.ok))
                stop("objective function returned NaN")
        }
        stopifnot(is.numeric(val))
        stopifnot(length(val) == 1)
        as.numeric(val)
    }
    if (is.list(objective))
        lapply(objective, objFunVal1)
    else
        objFunVal1(objective)
}

objFunVal.tf <-
objFunVal.hydromad <-
    function(x, objective = hydromad.getOption("objective"),
             ..., all = FALSE, nan.ok = FALSE)
{
    model <- x
    ## these can be referred to in `objective`
    X <- fitted(x, all = all)
    if (length(X) == 0)
        stop("fitted() returned nothing")
    delayedAssign("Q", observed(x, all = all))
    delayedAssign("U", fitted(x, all = all, U = TRUE))
    delayedAssign("P", observed(x, all = all, item = "P"))
    isValidModel <- isValidModel(x)
    objFunVal1 <- function(obj)
    {
        if (!isValidModel)
            return(NA_real_)
        if (inherits(obj, "formula"))
            obj <- obj[[2]]
        stopifnot(is.language(obj))
        val <- eval(obj)
        if (is.nan(val) && !nan.ok)
            stop("objective function returned NaN")
        stopifnot(is.numeric(val))
        stopifnot(length(val) == 1)
        as.numeric(val)
    }
    if (is.list(objective))
        lapply(objective, objFunVal1)
    else
        objFunVal1(objective)
}
