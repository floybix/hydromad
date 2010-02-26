## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


perfStats <-
    function(DATA, warmup = 0,
             which = c("rel.bias", "r.squared", "r.sq.sqrt", "r.sq.log", "r.sq.monthly"),
             na.action = na.pass)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    stopifnot(c("Q","X") %in% colnames(DATA))
    stopifnot(nrow(DATA) > warmup)
    warming <- seq_len(warmup)
    obs <- if (warmup>0) DATA[-warming,"Q"] else DATA[,"Q"]
    mod <- if (warmup>0) DATA[-warming,"X"] else DATA[,"X"]
    ans <- numeric()
    ## R Squared
    if ("r.squared" %in% which)
        ans[["r.squared"]] <- fitStat(obs, mod, p=2)
    ## R Squared square-root
    if ("r.sq.sqrt" %in% which) {
        ans[["r.sq.sqrt"]] <- fitStat(obs, mod, p=2, trans = sqrt)
    }
    ## R Squared log
    if ("r.sq.log" %in% which) {
        ans[["r.sq.log"]] <- fitStat(obs, mod, p=2, trans = log)
    }
    ## R Squared monthly
    if ("r.sq.monthly" %in% which) {
        ans[["r.sq.monthly"]] <- tsFitStat(obs, mod, p=2, aggr = 30)
    }
    ## Bias
    ## Rel. Bias
    if (any(c("bias", "rel.bias") %in% which)) {
        if ("bias" %in% which)
            ans[["bias"]] <- fitBias(obs, mod, rel = FALSE)
        if ("rel.bias" %in% which)
            ans[["rel.bias"]] <- fitBias(obs, mod)
    }
    ## U1, X1
    if (("U1" %in% which) && ("U" %in% colnames(DATA))) {
        U <- if (warmup>0) DATA[-warming,"U"] else DATA[,"U"]
        ans[["U1"]] <- cor(obs - mod, shiftWindow(U, -1), use="complete")
    }
    if ("X1" %in% which) {
        ans[["X1"]] <- cor(obs - mod, shiftWindow(mod, -1), use="complete")
    }
    ## Observed Gain
    if (("obsgain" %in% which) && ("U" %in% colnames(DATA)))
        ans[["obsgain"]] <- sum(DATA[-1,"Q"], na.rm=TRUE) /
            sum(DATA[,"U"], na.rm=TRUE)
    ans
}

objFunVal <-
    function(model, objective = hydromad.getOption("objective"), nan.ok = FALSE)
{
    if (inherits(objective, "formula"))
        objective <- objective[[2]]
    stopifnot(is.language(objective) || is.function(objective))
    if (!isValidModel(model))
        return(as.numeric(NA))
    ## these can be referred to in `objective`
    delayedAssign("Q", observed(model))
    delayedAssign("X", fitted(model))
    delayedAssign("U", fitted(model, U = TRUE))
    tmp <- eval(objective)
    if (nan.ok || is.nan(tmp))
        stop("objective function returned NaN")
    as.numeric(tmp)
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


zooFitStat <-
    function(obs, mod,
             ...)
{

}

tsFitStat <-
    function(obs, mod,
             ref = c("mean", "linear", "seasonal",
                      "blocks30", "blocks90", "blocks365"),
             trans = NULL,
             aggr = NULL,
             events = NULL,
             offset = identical(trans, log),
             na.action = na.pass,
             p = 2)
{
    obs <- as.ts(obs)
    mod <- as.ts(mod)
    stopifnot(length(obs) == length(mod))
    if (is.numeric(ref) && (length(ref) > 1))
        ref <- as.ts(ref)
    if (is.character(ref)) {
        if (length(ref) > 1)
            ref <- ref[1]
        if (any(grep("blocks", ref))) {
            block <- as.numeric(gsub("blocks", "", ref))
            ref <- "blocks"
        }
        ref <-
            switch(ref,
                   mean = mean(obs, na.rm = TRUE),
                   linear = {
                       ts(predict(lm(obs ~ time(obs)), newdata = time(obs)),
                          start = start(obs), frequency = frequency(obs))
                   },
                   seasonal = {
                       ## TODO
                       stop("'seasonal' not yet implemented")
                   },
                   blocks = {
                       levs <- rep(seq(1, ceiling(length(obs)/block)),
                                   each = block, length = length(obs))
                       vals <- tapply(obs, levs, FUN = mean, na.rm = TRUE)
                       ts(rep(vals, each = block, length = length(obs)),
                          start = start(obs), frequency = frequency(obs))
                   })
    }
    #if (length(ref) == 1) {
        ## if reference model is a single number (typically the mean)
        ## then we can do a quick version if not aggregating
    ref <- ts(rep(ref, length = length(obs)),
              start = start(obs), frequency = frequency(obs))
    stopifnot(length(ref) == length(obs))
    ## merge time series
    dat <- ts.intersect(obs = obs, mod = mod, ref = ref)
    if (NROW(dat) <= 1) {
        warning("merged time series have no data; incompatible times?")
        return(NA)
    }
    dat <- na.action(dat)
    if (NROW(dat) <= 1) {
        warning("time series have no data after 'na.action'")
        return(NA)
    }
    ## aggregation
    if (!is.null(aggr) && !is.null(events))
        stop("give at most one of 'aggr' and 'events'")
    if (!is.null(aggr)) {
        if (is.numeric(aggr))
            aggr <- list(ndeltat = aggr)
        if (!is.list(aggr))
            stop("unrecognised value of 'aggr'")
        dat <- do.call("aggregate", c(list(dat), aggr))
    }
    if (!is.null(events)) {
        if (!is.list(events))
            stop("unrecognised value of 'events'")
        FUN <- events$FUN
        if (is.null(FUN)) FUN <- "sum"
        events$FUN <- NULL
        ## compute events using first 2 series only (obs & mod)
        ev <- do.call("eventseq",
                      c(list(dat[,1:2]), events))
        ## apply FUN to events 'ev', in each series
        dat <- eventapply(dat, ev, FUN = FUN)
    }
    fitStat(dat[,"obs"], dat[,"mod"], ref = dat[,"ref"],
            trans = trans, offset = offset, p = p)
}


fitStat <-
    function(obs, mod,
             ref = NULL,
             trans = NULL,
             offset = identical(trans, log),
             p = 2)
{
    stopifnot(length(obs) == length(mod))
    if (!is.null(trans)) {
        ## offset = TRUE takes the observed 10%ile of non-zero values
        if (isTRUE(offset))
            offset <- quantile(obs[obs > 0], p=0.1, na.rm=TRUE)
        trans <- asSimpleFunction(trans, offset = offset)
        obs <- trans(obs)
        mod <- trans(mod)
        if (!is.null(ref))
            ref <- trans(ref)
    }
    ## default ref is the mean of transformed 'obs'
    if (is.null(ref))
        ref <- mean(obs, na.rm = TRUE)
    stopifnot(length(obs) == length(mod))
    ## calculate absolute error for model and reference
    ## and apply power p
    err <- abs(obs - mod) ^ p
    referr <- abs(obs - ref) ^ p
    ## only use pairwise common data
    ok <- complete.cases(err, referr)
    ans <- sum(err[ok]) / sum(referr[ok])
    1 - ans
}

fitBias <-
    function(obs, mod, relative.bias = TRUE, na.rm = TRUE)
{
    stopifnot(length(obs) == length(mod))
    err <- mod - obs
    ans <- mean(err, na.rm = na.rm)
    if (relative.bias)
        ans <- ans / mean(obs[!is.na(err)])
    ans
}

asSimpleFunction <- function(obj, offset = 0)
{
    if (is.character(obj))
        obj <- get(obj, mode = "function")
    if (inherits(obj, "formula")) {
        env <- environment(obj)
        funBody <- obj[[2]]
        obj <- function(x) NULL
        body(obj) <- funBody
        environment(obj) <- env
    }
    if (is.null(obj))
        obj <- force
    stopifnot(is.function(obj))
    if (offset != 0)
        return(function(x) obj(x + offset))
    return(obj)
}
