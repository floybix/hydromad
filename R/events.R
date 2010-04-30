##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

findThresh <-
    function(x, n, within = n %/% 20,
             mingap = 1, mindur = 1,
             below = FALSE, all = FALSE,
             trace = FALSE, optimize.tol = 0.1)
{
    stopifnot(n > 0)
    stopifnot(n * mindur + n * mingap <= NROW(x))
    x <- coredata(x)
    if (below) {
        x <- -x
    }
    ## return (difference from 'n' of) number of events for 'thresh'
    nDiffForThresh <- function(thresh) {
        ev <- eventseq(x, thresh = thresh, mingap = mingap,
                       mindur = mindur, all = all)
        newn <- nlevels(ev)
        if (trace)
            message(sprintf("thresh = %.3f, n = %d", thresh, newn))
        abs(newn - n)
    }

    ## do a quick run with a rough guess for the range
    rng <- quantile(x, (1 - (n * mindur / NROW(x))) ^ c(1, 3),
                     na.rm = TRUE, names = FALSE)
    if (trace)
        message("initial range guess: ", toString(signif(rng, 4)))
    res <- optimize(nDiffForThresh, interval = rng,
                    tol = optimize.tol)
    if (res$objective > within) {
        if (trace)
            message("initial run off by ", res$objective)
        ## expand possible range; only exclude minimum and maximum values
        rng2 <- quantile(x, range(ppoints(length(x))),
                        na.rm = TRUE, names = FALSE)
        res2 <- optimize(nDiffForThresh, interval = rng2,
                         tol = optimize.tol)
        if (res2$objective < res$objective)
            res <- res2
    }
    if (res$objective > within)
        warning("actual number of events differs from target (", n, ") by ",
                res$objective)
    thresh <- res$minimum
    if (below)
        thresh <- -thresh
    return(thresh)
}

eventseq <-
    function(x, thresh = 0, mingap = 1, mindur = 1, n = NA,
             below = FALSE, all = FALSE, inter = NA)
{
    if (!is.na(inter)) {
        .Deprecated(msg = "The 'inter' argument is deprecated. Use 'mingap'.")
        if (missing(mingap))
            mingap <- inter
    }
    if (!is.na(n)) {
        if (!missing(thresh))
            warning("'thresh' argument ignored since 'n' given")
        thresh <- findThresh(x, n = n, mingap = mingap, mindur = mindur,
                             below = below, all = all)
    }
    if (below) {
        x <- -x
        thresh <- -thresh
    }
    ## assume NAs are below threshold
    x[is.na(x)] <- -Inf
    ## find runs above threshold
    ## (runs continue while any series/column is above thresh)
    uruns <- if (length(dim(x)) == 2)
        rle(rowSums(coredata(x) > thresh) > 0)
    else rle(coredata(x) > thresh)
    ## find drops (between runs) whose length is too short
    if (mingap > 1) {
        nondrop <- with(uruns, values == FALSE & lengths < mingap)
        if (any(nondrop)) {
            ## (leave initial period alone)
            nondrop[1] <- FALSE
            ## set short drops to be part of the surrounding cluster
            uruns$values[nondrop] <- TRUE
            uruns <- rle(inverse.rle(uruns)) ## could be faster?
        }
    }
    ## find events whose length is too short
    if (mindur > 1) {
        nonev <- with(uruns, values == TRUE & lengths < mindur)
        if (any(nonev)) {
            ## (leave initial period alone)
            nonev[1] <- FALSE
            ## delete short events
            uruns$values[nonev] <- FALSE
            uruns <- rle(inverse.rle(uruns))
        }
    }
    ## assign unique numbers to events
    runvals <- uruns$values
    if (all) {
        ## assign numbers to events and inter-events
        ids <- seq_along(runvals)
        uruns$values <- ids
    } else {
        ## assign numbers to events, leave rest as NA
        ids <- seq_len(sum(runvals))
        uruns$values[runvals == TRUE] <- ids
        uruns$values[runvals == FALSE] <- NA
    }
    ## expand back into time series
    ev <- inverse.rle(uruns)
    ## the following should be same as ev <- factor(ev, ordered = TRUE)
    attr(ev, "levels") <- as.character(ids)
    class(ev) <- c("ordered", "factor")
    ## return events as a zoo with "factor" coredata
    ans <- zoo(ev, index(x))
    ans
}

eventapply <-
    function(X, events,
             FUN = sum, ..., by.column = TRUE, simplify = TRUE,
             TIMING = c("start", "middle", "end"))
{
    FUN <- match.fun(FUN)
    TIMING <- match.arg(TIMING)
    Xnames <- colnames(X)
    if (!is.vector(events) && !is.factor(events)) {
        ## merge series together (typically zoo or ts objects)
        cX <- cbind(X, events)
        events <- cX[,ncol(cX)]
        X <- cX[,-ncol(cX)]
        colnames(X) <- Xnames
    }
    ## need to handle functions returning vectors as well as scalars
    if (length(dim(X)) == 2) {
        ## multiple series
        if (by.column) {
            ## apply to each column separately
            ans <-
                lapply(as.data.frame(X), function(x) {
                    tmp <-
                        sapply(split(x, coredata(events), drop = TRUE),
                               FUN, ..., simplify = simplify)
                    if (is.matrix(tmp)) t(tmp) else tmp
                })
            ans <- as.matrix(do.call("data.frame", ans))
        } else {
            ## pass sub-period of multivariate series to function
            ans <-
                sapply(split(seq_len(NROW(X)), coredata(events), drop = TRUE),
                       function(ii) FUN(X[ii,,drop=FALSE], ...),
                       simplify = simplify)
            if (is.matrix(ans))
                ans <- t(ans)
        }

    } else {
        ## only one series
        ans <- sapply(split(X, coredata(events), drop = TRUE),
                      FUN, ..., simplify = simplify)
        if (is.matrix(ans))
            ans <- t(ans)
    }
    ## select time corresponding to each event
    timeIdxFn <- switch(TIMING,
                        start = function(x) x[1],
                        middle = function(x) x[ceiling(length(x)/2)],
                        end = function(x) x[length(x)])
    ev.index <- unlist(lapply(split(seq_len(NROW(X)), coredata(events), drop = TRUE),
                       timeIdxFn))
    ev.times <- time(X)[ev.index]
    if (simplify && !is.list(ans)) {
        return(zoo(as.matrix(ans), ev.times))
    } else {
        names(ans) <- format(ev.times)
        return(ans)
    }
}

preInterEventDuration <- function(x)
{
    if (inherits(x, "zoo"))
        stopifnot(is.regular(x, strict = TRUE))
    interCounter <- cumsum(is.na(coredata(x)))
    vals <- tapply(interCounter, x, FUN = head, 1)
    c(vals[1], diff(vals))
}

eventinfo <-
    function(X, events,
             FUN = mean, ...)
{
    stopifnot(inherits(X, "zoo"))
    ## TODO: allow FUN to return multiple values?
    xValues <- eventapply(X, events = events, FUN = FUN, ...)
    xLengths <- eventapply(X, events = events, FUN = length,
                           TIMING = "middle")
    midTimeComponents <- as.POSIXlt(time(xLengths))
    data.frame(Time = time(xValues),
               Month = midTimeComponents$mon + 1,
               Year = midTimeComponents$year + 1900,
               Value = coredata(xValues),
               Duration = coredata(xLengths),
               PreDuration = preInterEventDuration(events))
}
