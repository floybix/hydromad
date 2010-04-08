##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

event.clusters <- function(...) {
    .Deprecated("eventseq")
    eventseq(...)
}

eventAttributes <- function(...) {
    .Deprecated("eventinfo")
    eventinfo(...)
}


## based on the arguments of evd::clusters
eventseq <-
    function(x, thresh = 0, inter = 1, mindur = 1, below = FALSE)
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
    ids <- seq_len(sum(uruns$values))
    uruns$values[uruns$values == TRUE] <- ids
    uruns$values[uruns$values == FALSE] <- NA
    ## expand into long format
    ev <- inverse.rle(uruns)
    ## return events in the format of 'x', typically zoo or ts
    mostattributes(ev) <- attributes(x)
    attr(ev, "levels") <- ids
    ev
}

eventapply <-
    function(X, events = eventseq(X, thresh, inter, mindur, below),
             FUN = sum, ..., by.column = TRUE, simplify = TRUE,
             TIMING = c("start", "middle", "end"),
             thresh = 0, inter = 1, mindur = 1, below = FALSE)
{
    TIMING <- match.arg(TIMING)
    force(events)
    Xnames <- colnames(X)
    if (!is.vector(events) && !is.factor(events)) {
        ## merge series together (typically zoo or ts objects)
        cX <- cbind(X, events)
        events <- cX[,ncol(cX)]
        X <- cX[,-ncol(cX)]
        colnames(X) <- Xnames
    }
    ## extract only single numeric items if going to simplify
    .FUN <- FUN
    if (simplify)
        .FUN <- function(...) {
            tmp <- FUN(...)
            tmp <- tmp[unlist(lapply(tmp, function(z) {
                is.numeric(z) && !is.matrix(z) &&
                (length(z) == 1)
            }))]
            unlist(tmp)
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
                               .FUN, ..., simplify = simplify)
                    if (is.matrix(tmp)) t(tmp) else tmp
                })
            ans <- as.matrix(do.call("data.frame", ans))
        } else {
            ## pass sub-period of multivariate series to function
            ans <-
                sapply(split(seq_len(NROW(X)), coredata(events), drop = TRUE),
                       function(ii) .FUN(X[ii,,drop=FALSE], ...),
                       simplify = simplify)
            if (is.matrix(ans))
                ans <- t(ans)
        }

    } else {
        ## only one series
        ans <- sapply(split(X, coredata(events), drop = TRUE),
                      .FUN, ..., simplify = simplify)
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
    function(X, events = eventseq(X, thresh, inter, mindur, below),
             FUN = mean, ...,
             thresh = 0, inter = 1, mindur = 1, below = FALSE)
{
    stopifnot(inherits(X, "zoo"))
    force(events)
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
