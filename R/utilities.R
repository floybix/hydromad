## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## clean up regular zoo time series
tsCleanup <-
    function(x,
             na.action=na.pass,
             neg.rm=FALSE,
             return.info=TRUE,
             na.verb="removed")
{
    stopifnot(inherits(x, "zooreg"))
    ## trim off extraneous values at ends
    n <- NROW(x)
    x <- na.trim(x)
    trimmed <- n - NROW(x)
    ## detect and remove duplicated dates
    dup <- duplicated(index(x))
    ndups <- length(unique(index(x)[dup])) #sum(duplicated(index(x)))
    if (any(dup)) x <- aggregate(x, force, tail, 1)
    ## are any dates missing?
    regular <- is.regular(x, strict=TRUE)
    missdates <- 0
    if (!regular) {
        n <- NROW(x)
        xts <- as.ts(x)
        x <- zooreg(coredata(xts), start=start(x), #end=end(x),
                    frequency=frequency(x))
        #x <- byDays(as.ts(x))
        missdates <- NROW(x) - n
    }
    special <- is.infinite(x) | is.nan(x)
    nspecials <- sum(special)
    nnegs <- NA
    if (neg.rm) {
        neg <- (x < 0)
        nnegs <- sum(neg, na.rm=TRUE)
        if (nnegs > 0) x[neg] <- NA
    }
    missing <- sum(!complete.cases(x))
    missingPct <- round(100 * missing / NROW(x), digits=1)
    if (missingPct >= 1) missingPct <- round(missingPct)
    x <- na.action(x)
    finalmissing <- sum(!complete.cases(x))
    na.act <- missing - finalmissing
    info <- as.list(data.frame(missing, missingPct, finalmissing, na.act,
                                  trimmed, ndups, regular, missdates,
                                  nspecials, nnegs))
    desc <- paste(sep="", collapse="",
                  c(
                    if (neg.rm && (nnegs > 0)) c(nnegs, " negative values were removed. "),
                    if (missing == 0) "No time steps had missing values. "
                    else c(missing,
                           " time steps ",
                           " (", missingPct, " percent)",
                           " had missing values",
                           if (missdates > 0) c(" (", missdates, " of these were missing times)"),
                           if (na.act > 0) c(
                                             if (finalmissing == 0) " and all were " else
                                             c(", of which ", na.act, " were "),
                                             na.verb
                                             ),
                           ". "),
                    if (ndups > 0) c(ndups, " times were duplicated ",
                                     "(only the last value of each was kept). "),
                    if (nspecials > 0) c("Note: ", nspecials, " values were Inf or NaN.")
                    ))
    if (return.info) {
        attr(x, "cleanup.info") <- info
        attr(x, "cleanup.desc") <- desc
    }
    #return(list(object=x, info=info, description=desc))
    x
}


autocorrTime <-
    function(x, order = c(2,1), fallTo = 0.001)
{
    stopifnot(length(order) == 2)
    stopifnot((0 < fallTo) && (fallTo < 1))
    mod <- suppressWarnings(arima(x, order = c(order[1], 0, order[2]),
                                  include.mean = FALSE))
    pars <- coef(mod)
    ar <- pars[seq_len(order[1])]
    ma <- pars[seq_len(order[2]) + order[1] ]
    ans <- match(TRUE, ARMAtoMA(ar, ma, lag.max = NROW(x)/2) < fallTo)
    if (is.na(ans))
        ans <- NROW(x) / 2
    ans <- max(5, ans)
    ans
}

lagbind <-
    function(x, lags, all=TRUE, dframe=FALSE)
{
    bind.fn <- if (all) ts.union else ts.intersect
    names(lags) <- paste("t", ifelse(lags<0, "-", "+"),
                         abs(lags), sep="")
    do.call(bind.fn, c(lapply(lags, function(k) lag(x, k)),
                       list(dframe=dframe)))
}

shiftWindow <-
    function(x, delay, and.lag = FALSE, fill = NA)
{
    ## this is faster than a full ts.intersect - just shift by delay time steps
    ## a positive delay shifts the window back in time,
    ## a negative delay shifts the window forward in time
    if (delay == 0) return(x)
    attribs <- attributes(x)
    if (NCOL(x) == 1) {
        res <- if (delay > 0)
            c(rep(fill, delay), x[-(length(x)-seq(delay)+1)])
        else c(x[-seq(abs(delay))], rep(fill, abs(delay)))
    } else {
        res <- if (delay > 0)
            x[ c(rep(fill, delay), seq(1, NCOL(x) - delay)), ]
        else x[ c(seq(abs(delay)+1, NCOL(x)), rep(fill, abs(delay))), ]
    }
    attributes(res) <- attribs
    if (and.lag == FALSE) {
        ## already implicitly lagged above, the following keeps values in sync:
        res <- lag(res, delay)
    }
    res
}

stripWarmup <-
    function(x, warmup)
{
    if (length(x) == 0) return(x)
    if (warmup == 0) return(x)
    stopifnot(warmup >= 0)
    stopifnot(warmup < length(x))
    if (!is.ts(x)) x <- as.ts(x)
    freq <- frequency(x)
    t_warm <- tsp(x)[1] + (warmup) / freq
    t_end <- tsp(x)[2]
    ## following equivalent to return(window(x, start=t_warm)) but faster(?)
    if (NCOL(x) == 1)
        x <- x[-(1:warmup)]
    else 	x <- x[-(1:warmup),]
    ts(x, t_warm, t_end, freq)
}

byDays <- function(x)
{
    x <- as.zoo(x)
    index(x) <- as.Date(index(x), origin="1970-01-01")
    x
}
byYears <- function(x)
{
    x <- as.zoo(x)
    index(x) <- as.yearmon(index(x))
    x
}
bySecs <- function(x, tz="")
{
    x <- as.zoo(x)
    origin <- ISOdate(1970,1,1,0,0,0)
    attr(origin, "tzone") <- tz
    index(x) <- origin + as.numeric(index(x))
    x
}

observed <- function(object, ...)
    UseMethod("observed")

observed.default <- function(object, ...)
    fitted(object, ...) + residuals(object, ...)


numericSummary <- function(x) {
    x <- as.data.frame(x)
    foo <- sapply(x, function(z) {
        ans <- c(fivenum(z, na.rm=TRUE),
                 mean(z, na.rm=TRUE))
        names(ans) <- c("Min.", "1st Qu.", "Median",
                        "3rd Qu.", "Max.", "Mean")
        ## place mean after median
        ans <- ans[c(1:3,6,4:5)]
        nas <- sum(is.na(z))
        ans["NA's"] <- nas
        ans
    })
    #if (all(foo["NA's",] == 0)) ...
    foo
}



## TODO remove this when latticeExtra 0.6-7 released
simpleSmoothTs <- 
    function(x, width = NROW(x) %/% 10 + 1,
             c = 1, sides = 2, circular = FALSE,
             kern = kernel("daniell", rep(floor((width/sides)/sqrt(c)), c)))
{
    if (sides == 2) {
        ii <- -kern$m:kern$m
        filter <- kern[ii]
    } else if (sides == 1) {
        ii <- -kern$m:0
        filter <- kern[ii] / sum(kern[ii]) ## normalise
    } else stop("unrecognised value of 'sides'")
    xf <- x
    xf[] <- filter(as.matrix(x), filter, sides = sides, circular = circular)
    xf
}
