## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

tsFitStat <-
    function(obs, mod, ref = NULL, ...,
             na.action = na.pass,
             aggr = NULL, events = NULL)
{
    if (is.null(ref)) {
        ref <- mean(obs, na.rm = TRUE)
    }
    if (length(ref) == 1) {
        ## if reference model is a single number (typically the mean)
        ## turn it into a time series like 'obs'
        tmp <- ref
        ref <- obs
        ref[] <- tmp
    }
    ## merge time series
    dat <- cbind(obs = obs, mod = mod, ref = ref)
    if (NROW(dat) <= 1) {
        warning("merged time series have no data; incompatible times?")
        return(NA_real_)
    }
    dat <- na.action(dat)
    if (NROW(dat) <= 1) {
        warning("time series have no data after 'na.action'")
        return(NA_real_)
    }
    ## aggregation
    if (!is.null(aggr) && !is.null(events))
        stop("give at most one of 'aggr' and 'events'")
    if (!is.null(aggr)) {
        if (is.numeric(aggr))
            aggr <- list(ndeltat = aggr)
        if (is.character(aggr)) {
            dat <- aggregate(dat, by = cut(time(dat), aggr))
        } else {
            if (!is.list(aggr))
                stop("unrecognised value of 'aggr'")
            dat <- do.call("aggregate", c(alist(dat), aggr))
        }
    }
    if (!is.null(events)) {
        if (!is.list(events))
            stop("unrecognised value of 'events'")
        FUN <- events$FUN
        if (is.null(FUN)) FUN <- "sum"
        events$FUN <- NULL
        ## compute events using first 2 series only (obs & mod)
        ev <- do.call("eventseq",
                      c(alist(dat[,"obs"]), events))
        ## apply FUN to events 'ev', in each series
        dat <- eventapply(dat, ev, FUN = FUN)
    }
    fitStat(dat[,"obs"], dat[,"mod"], ref = dat[,"ref"], ...)
}


fitStat <-
    function(obs, mod, ref = NULL, p = 2,
             trans = NULL, offset = identical(trans, log),
             negatives.ok = FALSE)
{
    if (!identical(attributes(obs), attributes(mod))) {
        warning("attributes of 'obs' and 'mod' are not identical; need to use tsFitStat?")
    }
    obs <- coredata(obs)
    mod <- coredata(mod)
    stopifnot(length(obs) == length(mod))
    ## only use pairwise common data
    ok <- complete.cases(obs, mod)
    if (length(ref) > 1) {
        ref <- coredata(ref)
        ok <- ok & !is.na(ref)
        ref <- ref[ok]
    }
    obs <- obs[ok]
    mod <- mod[ok]
    if (length(obs) == 0)
        return(NA_real_)
    ## negative values can cause errors with log()
    if (negatives.ok == FALSE) {
        obs <- pmax(obs, 0)
        mod <- pmax(mod, 0)
        if (!is.null(ref))
            ref <- pmax(ref, 0)
    }
    if (!is.null(trans)) {
        ## offset = TRUE takes the observed 10%ile of non-zero values
        if (isTRUE(offset)) {
            offset <- quantile(obs[obs > 0], p = 0.1)
            if (!is.finite(offset)) offset <- 0
        }
        trans <- asSimpleFunction(trans, offset = offset)
        obs <- trans(obs)
        mod <- trans(mod)
        if (!is.null(ref))
            ref <- trans(ref)
    }
    ## check again, just in case 'trans' changed the length
    stopifnot(length(obs) == length(mod))
    ## default ref is the mean of transformed 'obs'
    if (is.null(ref)) {
        ref <- mean(obs, na.rm = TRUE)
    }
    ## calculate absolute error for model and reference
    ## and apply power p
    err <- abs(obs - mod) ^ p
    referr <- abs(obs - ref) ^ p
    ans <- sum(err) / sum(referr)
    1 - ans
}

fitBias <-
    function(obs, mod, relative.bias = TRUE)
{
    ## will signal a warning if not same length (or a multiple):
    err <- mod - obs
    ok <- !is.na(err)
    ans <- mean(err[ok])
    if (relative.bias)
        ans <- ans / mean(obs[ok])
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
