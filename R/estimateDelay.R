## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


estimateDelay <-
    function(DATA = data.frame(U=, Q=),
             rises = TRUE, rank = FALSE, n.estimates = 1,
             lag.max = hydromad.getOption("max.delay"),
             na.action = na.exclude, plot = FALSE, main = NULL,
             ...)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    ## which column represents streamflow?
    iQ <- 1
    if ("Q" %in% colnames(DATA)) {
        iQ <- which("Q" == colnames(DATA))[1]
    }
    Q <- DATA[,iQ]
    U <- DATA[,if (iQ==1) 2 else 1]
    ylab <- "CCF"
    do.rises <- rises
    if (do.rises) {
        ## backwards difference, i.e. rises are placed at their end time
        rises <- function(x) { x[] <- c(0, pmax(diff(x), 0)); x }
        Q <- rises(Q)
    }
    if (rank) {
        Q[] <- rank(zapsmall(Q), na.last="keep")
        U[] <- rank(zapsmall(U), na.last="keep")
        ylab <- "rank CCF"
    }
    if (is.null(main)) {
        main <- "Cross-correlation"
        if (do.rises) main <- paste(main, "with rises only")
    }
    if (sum(complete.cases(Q, U)) == 0)
        return(NA)
    ans <- ccf(Q, U, lag.max = lag.max, na.action = na.action,
               plot = plot, main = main, ...)
    ##ans$lag[which.max(ans$acf)]
    est <- ans$lag[order(ans$acf, decreasing=TRUE)]
    #if (est[1] < 0) warning("delay estimate is negative, check your data")
    ## omit any negative delays or NAs
    invalid <- is.na(est) | (est < 0)
    #invalid[1] <- FALSE
    est <- est[!invalid]
    est <- est[seq(n.estimates)]
    ## convert from "basic units" to time steps
    est <- est * frequency(Q)
    est
}

