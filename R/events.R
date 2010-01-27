## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##

## based on the arguments of evd::clusters
event.clusters <- function(x, thresh = 0, inter = 1, mindur = 1, below = FALSE)
                                        #, ulow = -Inf, rlow = 1)
{
    if (below) {
        x <- -x
        thresh <- -thresh
    }
    ## assume NAs are below threshold
    x[is.na(x)] <- -Inf
    ## find runs above threshold
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
    ## find clusters whose length is too short
    if (mindur > 1) {
        nonclus <- with(uruns, values == TRUE & lengths < mindur)
        if (any(nonclus)) {
            ## (leave initial period alone)
            nonclus[1] <- FALSE
            ## delete short clusters
            uruns$values[nonclus] <- FALSE
            uruns <- rle(inverse.rle(uruns))
        }
    }
    ## assign unique numbers to clusters
    uruns$values[uruns$values == TRUE] <- seq_len(sum(uruns$values))
    uruns$values[uruns$values == FALSE] <- NA
    ## return clusters as a vector or ts, to be used in tapply etc
    clust <- inverse.rle(uruns)
    mostattributes(clust) <- attributes(x)
    clust
}

eventapply <-
    function(X,
             clusters = event.clusters(X, thresh=thresh, inter=inter, mindur=mindur),
             FUN = sum, ...,
             thresh = 0, inter = 1, mindur = 1,
             TIMING = c("start", "middle", "end"),
             by.column = TRUE)
{
    TIMING <- match.arg(TIMING)
    force(clusters)
    ## handle functions returning vectors
    if (length(dim(X)) == 2) {
        ## multiple series
        ## TODO: test this
        if (by.column) {
            ## apply to each column separately
            ans <- apply(X, 2, function(x) {
                tmp <- tapply(x, clusters, FUN, ...)
                if (is.list(tmp))
                    tmp <- do.call(rbind, tmp)
                tmp
            })
        } else {
            ## pass sub-period of multivariate series to function
            ans <- tapply(seq(NROW(X)), clusters,
                          function(ii) FUN(X[ii,,drop=FALSE], ...))
            if (is.list(ans))
                ans <- do.call(rbind, ans)
        }

    } else {
        ## only one series
        ans <- tapply(X, clusters, FUN, ...)
        if (is.list(ans))
            ans <- do.call(rbind, ans)
    }
    ## select time corresponding to each cluster
    timeIdxFn <- switch(TIMING,
                        start = function(x) x[1],
                        middle = function(x) x[ceiling(length(x)/2)],
                        end = function(x) x[length(x)])
    clust.index <- tapply(seq(NROW(X)), clusters, timeIdxFn)
    zoo(as.matrix(ans), time(X)[clust.index])
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
             clusters = event.clusters(X, thresh=thresh, inter=inter, mindur=mindur),
             FUN = mean, ...,
             thresh = 0, inter = 1, mindur = 1)
{
    stopifnot(inherits(X, "zoo"))
    force(clusters)
    xValues <- eventapply(X, clusters = clusters, FUN = FUN, ...)
    xLengths <- eventapply(X, clusters = clusters, FUN = length,
                           TIMING = "middle")
    data.frame(Date = time(xValues),
               Month = as.POSIXlt(time(xLengths))$mon + 1,
               Value = coredata(xValues),
               Duration = coredata(xLengths),
               DryPeriod = preInterEventDuration(clusters))
}
