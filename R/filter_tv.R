## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


## time-varying recursive filter
filter_tv <-
    function(x, a, init = 0)
{
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(init))
    stopifnot(length(x) == length(a))
    ## skip over missing values (maintaining the state y[i-1])
    bad <- is.na(x) | is.na(a)
    x[bad] <- 0
    a[bad] <- 1
    if (ihacres.getOption("pure.R.code")) {
        y <- x
        y[1] <- x[1] + a[1] * init
        for (i in 2:length(x)) {
            y[i] <- x[i] + a[i] * y[i-1]
        }
    } else {
        y <- .C(ar1_tv,
                as.double(x),
                as.double(a),
                as.integer(length(x)),
                as.double(init),
                out = double(length(x)),
                DUP=FALSE, PACKAGE="ihacreslab")$out
        ## make it a time series object again
        attributes(y) <- attributes(x)
    }
    ## re-insert missing values
    y[bad] <- NA
    y
}
