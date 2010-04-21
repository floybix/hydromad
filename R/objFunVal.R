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
         "r.sq.log" = ~ fitStat(Q, X, trans = log),
         "r.sq.diff" = ~ fitStat(diff(Q), diff(X)),
         "r.sq.rank" =
         ~ fitStat(Q, X, trans = function(x)
                   rank(zapsmall(x, digits = 3), ties = "min", na.last = "keep")),
         "r.sq.monthly" =
         ~ tsFitStat(Q, X, aggr = list(by = cut, breaks = "months", FUN = sum)),
         "persistence" = ~ tsFitStat(Q, X, ref = lag(Q, -1)),
         "r.sq.seasonal" =
         ~ tsFitStat(Q, X, ref = ave(Q, months(time(Q)))),
         "r.sq.vs.P" =
         ~ fitStat(Q, X, ref = { rx <- filter(P, ar(Q, demean=FALSE)$ar, "r");
                                 rx * mean(Q, na.rm = TRUE) /
                                     mean(rx, na.rm = TRUE) }),
         "events.med5sum" =
         ~ tsFitStat(Q, X, events = list(thresh = median(Q, na.rm = TRUE), inter = 5, FUN = sum)),
         "events.med5max" =
         ~ tsFitStat(Q, X, events = list(thresh = median(Q, na.rm = TRUE), inter = 5, FUN = max)),
         "events.med5min" =
         ~ tsFitStat(Q, X, events = list(thresh = median(Q, na.rm = TRUE), below = TRUE, mindur = 5, FUN = min)),
         "abs.err" = ~ mean(abs(X - Q), na.rm = TRUE),
         "RMSE" = ~ sqrt(mean((Q - X)^2, na.rm = TRUE)),
         "U1" = ~ cor(Q - X, shiftWindow(U, -1), use="complete"),
         "X1" = ~ cor(Q - X, shiftWindow(X, -1), use="complete")
         )
}

objFunVal <- function(x, objective, ...)
    UseMethod("objFunVal")

objFunVal.default <-
    function(x, objective = hydromad.getOption("objective"),
             ..., nan.ok = FALSE)
{
    ## these can be referred to in `objective`
    delayedAssign("Q", x[,"Q"])
    delayedAssign("X", x[,"X"])
    delayedAssign("U", x[,"U"])
    delayedAssign("P", x[,"P"])
    objFunVal1 <- function(obj)
    {
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

objFunVal.hydromad <-
    function(x, objective = hydromad.getOption("objective"),
             ..., all = FALSE, nan.ok = FALSE)
{
    ## these can be referred to in `objective`
    model <- x
    delayedAssign("Q", observed(x, all = all))
    delayedAssign("X", fitted(x, all = all))
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

bestByObjFun <-
    function(models,
             objective = hydromad.getOption("objective"),
             maximum = TRUE,
             orStop = FALSE)
{
    stopifnot(is.list(models) && (length(models) > 0))
    ## screen out invalid models
    ok <- sapply(models, isValidModel)
    if (!any(ok)) {
        if (orStop)
            stop("found no valid models")
        ## otherwise just return an invalid model
        return(models[[1]])
    }
    models <- models[ok]
    if (length(models) == 1)
        return(models[[1]])
    objVals <-
        lapply(models, objFunVal, objective = objective)
    objVals <- unlist(objVals) * ifelse(maximum, 1, -1)
    models[[ which.max(objVals) ]]
}
