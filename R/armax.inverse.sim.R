## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


armax.inverse.sim <-
    function(Q, P = NULL,
             pars = c(a_1 = 0, b_0 = 1),
             delay = 0, Xs_0 = 0, Xq_0 = 0,
             rain.factor = 1.1,
             rises.only = FALSE,
             use.Qm = TRUE,
             use.fft.method = FALSE,
             constrain.fft = TRUE,
             mass.balance = use.fft.method,
             scale.window = NA)
{
    if (!is.ts(Q)) Q <- as.ts(Q)
    pars <- tfParsConvert(pars, "a,b")
    tfParsCheck(pars)
    ## extract first parts of parameter names (before underscore)
    parSymbols <- gsub("_.*", "", names(pars))
    a <- pars[parSymbols == "a"]
    b <- pars[parSymbols == "b"]
    n <- length(a)
    m <- length(b) - 1
                                        #if (constrain.zeros) {
                                        #    if (!is.null(P)) {
                                        #        zeros <- mean(zapsmall(P, 2) == 0, na.rm = TRUE)
                                        #    } else {
                                        #        zeros <- mean(zapsmall(diff(Q), 2) <= 0, na.rm = TRUE)
                                        #    }
                                        #}
    ## set NAs in flow to zero (insert them again afterwards)
    isna <- is.na(Q)
    Q[isna] <- 0
    ## initialise P
    if (!is.null(P)) {
        if (!is.ts(P)) P <- as.ts(P)
        P <- P * rain.factor
        ## apply delay
        if (delay > 0) {
            P <- shiftWindow(P, delay = delay, and.lag = TRUE)
        }
        P <- window(P, start = start(Q), end = end(Q), extend = TRUE)
    } else {
        P <- Q
        P[] <- 1e12                     #Inf
    }
    if (rises.only) {
        ## constrain U to time steps with rising flow
        diffrising <- (diff(Q) > 0)
        ## backwards difference, i.e. rises are placed at their end time
        rising <- c(FALSE, diffrising)
        P[!rising] <- 0
    }
    ## restrict flow on missing rain days
    P[is.na(P)] <- 0
    stopifnot(length(P) == length(Q))
    ## initialise U
    U <- Q
    U[] <- 0

    if (use.fft.method) {
        ## fourier deconvolution method
        bb <- b[-1] / b[1]
        maform <- ARMAtoMA(a, bb, lag.max = length(Q))
        U <- wiener(Q, maform, N = 0)
        U <- pmax(U, 0)
                                        #U[U <= quantile(U, zeros, na.rm = TRUE)] <- 0
        if (constrain.fft)
            U <- pmin(U, P)
    } else

    if (hydromad.getOption("pure.R.code")) {
        ## slow version in R for cross-checking
        Qm <- Q
        for (t in seq(max(n,m)+1, length(Q))) {
            Ut <- Q[t]
            if (n > 0)
                Ut <- Ut - sum(a * Qm[t - seq_len(n)])
            if (m > 0)
                Ut <- Ut - sum(b[-1] * U[t - seq_len(m)])
                                        #if (n >= 1) Ut <- Ut - a[1] * Qm[t-1]
                                        #if (n >= 2) Ut <- Ut - a[2] * Qm[t-2]
                                        #if (n >= 3) Ut <- Ut - a[3] * Qm[t-3]
                                        #if (m >= 1) Ut <- Ut - b[2] * U[t-1]
                                        #if (m >= 2) Ut <- Ut - b[3] * U[t-2]
                                        #if (m >= 3) Ut <- Ut - b[4] * U[t-3]
            Ut <- Ut / b[1]
            U[t] <- max(0, min(Ut, P[t]))
            if (use.Qm) {
                Qm[t] <- (sum(a * Qm[t - seq_len(n)]) +
                          sum(b * U[t - 0:m]))
                                        #if (n >= 1) Qmt <- Qmt + a[1] * Qm[t-1]
                                        #if (n >= 2) Qmt <- Qmt + a[2] * Qm[t-2]
                                        #if (n >= 3) Qmt <- Qmt + a[3] * Qm[t-3]
                                        #if (m >= 0) Qmt <- Qmt + b[1] * U[t]
                                        #if (m >= 1) Qmt <- Qmt + b[2] * U[t-1]
                                        #if (m >= 2) Qmt <- Qmt + b[3] * U[t-2]
                                        #if (m >= 3) Qmt <- Qmt + b[4] * U[t-3]
                                        #Qm[t] <- Qmt
            }
        }
    } else {
        U <- .C(inverse_filter,
                as.double(Q[]), ## force copy
                as.integer(length(Q)),
                as.double(a),
                as.double(b),
                as.integer(n),
                as.integer(m),
                as.integer(use.Qm),
                as.double(P),
                U = double(length(Q)),
                DUP=FALSE, PACKAGE="hydromad")$U
        ## make it a time series object again
        attributes(U) <- attributes(Q)
    }
    U[isna] <- NA
    ## scale to ensure mass balance with Q
    if (mass.balance) {
        if (is.na(scale.window)) {
            U <- U * sum(Q[!isna]) / sum(U[!isna])
        } else {
            sc <- simpleSmoothTs(cbind(Q, U), width = scale.window)
            sc <- sc[,"Q"] / sc[,"U"]
            sc <- na.locf(na.locf(sc, na.rm = FALSE), fromLast = TRUE, na.rm = FALSE)
            U <- U * sc
        }
    }
    ## inverse.sim should be in sync with P not Q,
    ## so delay still applies
    if (delay > 0) {
        U <- shiftWindow(U, delay = -delay, and.lag = TRUE)
    }
    U
}


## based on http://en.wikipedia.org/w/index.php?title=Wiener_filter&oldid=211541130
wiener <- function(y, h, N, gamma = 0)
{
    length(h) <- length(y)
    h[is.na(h)] <- 0
    y[is.na(y)] <- 0
    H <- fft(h)
    Y <- fft(y)
    S <- abs(Y)^2 / length(y)^2

    ## experimental threshold to reduce sensitivity
    ## based on http://cnx.org/content/m13144/latest/wienerFilter.m
    ## TODO: is this ever useful?
    zap <- Mod(H) < gamma
    if (any(zap))
        H[zap] <- gamma * Mod(H[zap]) / H[zap]

    G <- (Conj(H) * S) /
        (abs(H)^2 * S + N)

    x <- Re(fft(G * Y, inverse=TRUE) / length(y))
    x
}
