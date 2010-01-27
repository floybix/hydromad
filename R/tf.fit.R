## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


tf.fit <-
    function(..., method = ihacres.getOption("uh.method"))
{
    stopifnot(is.character(method))
    fnName <- paste("tf", method, "fit", sep = ".")
    uhcall <- quote(tf.fit(...))
    uhcall[[1]] <- as.symbol(fnName)
    eval(uhcall)
}

defaultPrefilters <- function()
{
    alpha_q <- c(0.01, 0.2)
    alpha_s <- c(0.9, 0.95, 0.98, 0.8)

    alpha <- expand.grid(alpha_q = alpha_q, alpha_s = alpha_s)

    a1 <- with(alpha, alpha_s + alpha_q)
    a2 <- with(alpha, -(alpha_s * alpha_q))
    aa <- data.frame(a_1 = a1, a_2 = a2)
    split(as.matrix(aa), 1:NROW(aa))
}

makePrefilter <-
    function(DATA = list(Q=),
             order = ihacres.getOption("order"),
             pureAR = FALSE,
             na.action = na.exclude)
{
    ## get data into the right form
    if (is.list(DATA)) DATA <- do.call(ts.intersect, lapply(DATA, as.ts))
    if (!is.ts(DATA)) DATA <- as.ts(DATA)
    if (NCOL(DATA) > 1)
        DATA <- DATA[,"Q"]
    DATA <- na.action(DATA)
    ## check values
    stopifnot(length(order) == 2)
    order[1] <- max(1, order[1])
    ## generate prefilter from data
    ## pure AR method is very fast
    if (pureAR)
        return(ar(DATA, order.max = order[1], demean = FALSE,
                  na.action = na.action)$ar)
    ## otherwise fit an ARMA model
    mod <- arima(DATA, order = c(order[1], 0, order[2]+1),
                 include.mean = FALSE)
    coef(mod)[1:(order[1])]
}


alpha.init <-
    function(n = 2)
{
    stopifnot(n %in% 1:2)

    alpha_q <- c(0.01, 0.2)
    alpha_s <- c(0.9, 0.95, 0.98, 0.8)

    if (n == 1)
        return(data.frame(alpha_s = alpha_q))

    alpha <- expand.grid(alpha_q = alpha_q, alpha_s = alpha_s)
    alpha <- subset(alpha, alpha_q <= alpha_s)
    return(alpha)
    #unlist(alpha[k,])
}

a.init <-
    function(n = 2)
{
    alpha <- alpha.init(n = n)
    if (n == 1) return(alpha)
    a1 <- with(alpha, alpha_s + alpha_q)
    a2 <- with(alpha, -(alpha_s * alpha_q))
    data.frame(a_1 = a1, a_2 = a2)
}



estimateOrder.heuristic <-
    function(DATA=list(U=, Q=), delay=NULL)
{
    ## get data into the right form
    if (is.list(DATA)) DATA <- do.call(ts.intersect, lapply(DATA, as.ts))
    if (!is.ts(DATA)) DATA <- as.ts(DATA)
    if (!("U" %in% colnames(DATA)) && ("P" %in% colnames(DATA)))
        colnames(DATA) <- gsub("P", "U", colnames(DATA))
    stopifnot(c("U","Q") %in% colnames(DATA))
    U <- DATA[,"U"]
    Q <- DATA[,"Q"]
    if (is.null(delay))
        delay <- estimateDelay(DATA, plot=FALSE, na.action=na.exclude)
    ## start with most basic model: single exponential
    order <- c(n=1, m=0)
    ##forcing <- which(U > 0) + delay
    ##isforcing <- seq_along(Q) %in% forcing

    ## use two exponential components?
    isforcing <- (U > 0)
    ## and up to delay steps after input
    for (d in seq_len(delay))
        isforcing <- isforcing | lag(isforcing, -1)
    ## re-align
    isforcing <- c(rep(NA, delay), isforcing)
    Q.falling <- Q
    Q.falling[isforcing] <- NA
    ## only consider another exponential if at least 25% non-zero flow
    zero.frac <- sum(Q.falling == 0, na.rm=TRUE) / length(Q.falling)
    if (zero.frac < 0.75) {
        rates <- Q.falling / lag(Q.falling, -1)
        rates <- rates[is.finite(rates)]
        rates <- rates[(0 < rates) & (rates < 1)]
        ##rrates <- recess[U == 0]
        ##recess <- ifelse(DATA[,"U"] == 0, DATA[,"Q"], NA)
        ##rrates <- recess / lag(recess, -1)
        rtimes <- -1 / log(rates)
        ##rtimesdev <- IQR(rtimes, na.rm=T)
        ##if (rtimesdev >= 4) ## inter-quartile range of decay time >= 4 timesteps

        ## if working with very little data, can not justify another exponential
        if (length(rtimes) > 100)
            if (quantile(rtimes, 0.5) > 2 * quantile(rtimes, 0.05))
                order <- order + 1
    }
    ## include instantaneous component?
    fall.0 <- which(diff(U > 0) == -1) + delay
    clearfalls <- (U[fall.0] == 0) & (U[fall.0+1] == 0) & (U[fall.0+2] == 0)
    fall.0 <- na.omit(fall.0[clearfalls])
    if (length(fall.0) < 30)
        return(order)
    rates.0 <- Q[fall.0+1] / Q[fall.0]
    rates.1 <- Q[fall.0+2] / Q[fall.0+1]
    ok <- is.finite(rates.0) & is.finite(rates.1)
    fall.0 <- fall.0[ok]
    rates.0 <- rates.0[ok]
    rates.1 <- rates.1[ok]
    rates.0.1 <- rates.1 / rates.0
    rate.ratio.m <- weighted.mean(rates.0.1, Q[fall.0])
    ##mrate.0 <- weighted.mean(rates.0, Q[fall.0])
    ##mrate.1 <- weighted.mean(rates.1, Q[fall.0])
    ##if (mrate.0 < 0.75 * mrate.1)
    if (rate.ratio.m < 0.75)
        order["m"] <- order["m"] + 1
    ##recover()
    order
}

estimateDelay <-
    function(DATA = data.frame(U=, Q=),
             rises=TRUE, rank=FALSE, n.estimates=1,
             lag.max = ihacres.getOption("max.delay"),
             na.action=na.exclude, plot=FALSE, main=NULL, ...)
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

