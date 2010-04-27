##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

eventseq <-
    function(x, thresh = 0, mingap = 1, mindur = 1,
             below = FALSE, all = FALSE, inter = NA)
{
    if (!is.na(inter)) {
        .Deprecated(msg = "The 'inter' argument is deprecated. Use 'mingap'.")
        if (missing(mingap))
            mingap <- inter
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
    ids <- seq_len(sum(runvals))
    uruns$values[runvals == TRUE] <- ids
    if (all) {
        ## assign negative numbers to inter-events
        inter.ids <- - seq_len(sum(runvals == FALSE))
        uruns$values[runvals == FALSE] <- inter.ids
        ids <- c(ids, inter.ids)
    } else {
        ## set inter-event periods to NA
        uruns$values[runvals == FALSE] <- NA
    }
    ## expand into long format
    ev <- inverse.rle(uruns)
    ## return events in the format of 'x', typically zoo or ts
    mostattributes(ev) <- attributes(x)
    class(ev) <- unique(c("eventseq", class(ev)))
    #attr(ev, "levels") <- ids
    ev
}

levels.eventseq <- function(x)
    unique(x[!is.na(coredata(x))])

eventapply <-
    function(X, events,
             FUN = sum, ..., by.column = TRUE, simplify = TRUE,
             TIMING = c("start", "middle", "end"))
{
    TIMING <- match.arg(TIMING)
    Xnames <- colnames(X)
    if (!is.vector(events) && !is.factor(events)) {
        ## merge series together (typically zoo or ts objects)
        cX <- cbind(X, events)
        events <- cX[,ncol(cX)]
        X <- cX[,-ncol(cX)]
        colnames(X) <- Xnames
    }
    ## extract only single numeric items if going to simplify
    ## TODO: need this?
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
