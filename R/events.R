##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## based on the arguments of evd::clusters
event.clusters <-
eventSeq <-
    function(x, thresh = 0, inter = 1, mindur = 1, below = FALSE)
                                        #, ulow = -Inf, rlow = 1)
{
    if (below) {
        x <- -x
        thresh <- -thresh
    }
    ## assume NAs are below threshold
    x[is.na(x)] <- -Inf
    ## find runs above threshold
    ## (runs are merged across columns)
    uruns <- if (length(dim(x)) == 2)
        rle(rowSums(coredata(x) > thresh) > 0)
    else rle(coredata(x) > thresh)
    ## find drops (between runs) whose length is too short
    if (inter > 1) {
        nondrop <- with(uruns, values == FALSE & lengths < inter)
        if (any(nondrop)) {
            ## (leave initial period alone)
            nondrop[1] <- FALSE
            ## set short drops to be part of the surrounding cluster
            uruns$values[nondrop] <- TRUE
            uruns <- rle(inverse.rle(uruns))
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
    uruns$values[uruns$values == TRUE] <- seq_len(sum(uruns$values))
    uruns$values[uruns$values == FALSE] <- NA
    ## expand into long format
    ev <- inverse.rle(uruns)
    ## return events in the format of 'x', typically zoo or ts
    mostattributes(ev) <- attributes(x)
    ev
}

eventapply <-
    function(X,
             events = eventSeq(X, thresh = thresh, inter = inter, mindur = mindur, below = below),
             FUN = sum, ...,
             thresh = 0, inter = 1, mindur = 1, below = FALSE,
             TIMING = c("start", "middle", "end"),
             by.column = TRUE)
{
    TIMING <- match.arg(TIMING)
    force(events)
    ## merge series together (typically zoo or ts objects)
    cX <- cbind(events, X)
    events <- cX[,1]
    X <- cX[,-1]
    ## handle functions returning vectors
    if (length(dim(X)) == 2) {
        ## multiple series
        ## TODO: test this
        if (by.column) {
            ## apply to each column separately
            ans <- apply(X, 2, function(x) {
                tmp <- tapply(x, events, FUN, ...)
                if (is.list(tmp))
                    tmp <- do.call(rbind, tmp)
                tmp
            })
        } else {
            ## pass sub-period of multivariate series to function
            ans <- tapply(seq(NROW(X)), events,
                          function(ii) FUN(X[ii,,drop=FALSE], ...))
            if (is.list(ans))
                ans <- do.call(rbind, ans)
        }

    } else {
        ## only one series
        ans <- tapply(X, events, FUN, ...)
        if (is.list(ans))
            ans <- do.call(rbind, ans)
    }
    ## select time corresponding to each event
    timeIdxFn <- switch(TIMING,
                        start = function(x) x[1],
                        middle = function(x) x[ceiling(length(x)/2)],
                        end = function(x) x[length(x)])
    ev.index <- tapply(seq(NROW(X)), events, timeIdxFn)
    zoo(as.matrix(ans), time(X)[ev.index])
}

preInterEventDuration <- function(x)
{
    if (inherits(x, "zoo"))
        stopifnot(is.regular(x, strict = TRUE))
    interCounter <- cumsum(is.na(coredata(x)))
    vals <- tapply(interCounter, x, FUN = head, 1)
    c(vals[1], diff(vals))
}

eventAttributes <-
    function(X,
             events = eventSeq(X, thresh = thresh, inter = inter, mindur = mindur, below = below),
             FUN = mean, ...,
             thresh = 0, inter = 1, mindur = 1, below = FALSE)
{
    stopifnot(inherits(X, "zoo"))
    force(events)
    xValues <- eventapply(X, events = events, FUN = FUN, ...)
    xLengths <- eventapply(X, events = events, FUN = length,
                           TIMING = "middle")
    data.frame(Date = time(xValues),
               Month = as.POSIXlt(time(xLengths))$mon + 1,
               Value = coredata(xValues),
               Duration = coredata(xLengths),
               DryPeriod = preInterEventDuration(events))
}
