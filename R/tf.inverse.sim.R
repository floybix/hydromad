## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


tf.fft.inverse.sim <-
    function(Q, ..., use.fft.method = TRUE)
{
    tf.inverse.sim(Q, ..., use.fft.method = use.fft.method)
}

tf.inverse.sim <-
    function(Q, P = NULL,
             pars = c(tau_s = 0, tau_q = 0, v_s = 1, v_3 = 0),
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
    if ("lambda" %in% parSymbols) {
        use.fft.method <- FALSE
    }
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
            P <- shiftWindow(P, delay = delay)
            ## (and fix up time series attributes)
            P <- lag(P, -delay)
        }
    } else {
        P <- Q
        P[] <- 1e12 #Inf
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

    ######TODO
    if (TRUE || hydromad.getOption("pure.R.code")) {
        ## slow version in R for cross-checking
        if ("lambda" %in% parSymbols) {
            pars.tv <- tfParsConvert(pars, "tau,v")
            pars.alpha <- tfParsConvert(pars, "alpha,beta")
            alpha_s <- pars.alpha["alpha_s"]
            alpha_q <- pars.alpha["alpha_q"]
            if (is.na(alpha_q)) alpha_q <- 0
            v_s_0 <- pars.tv["v_s"]
            lambda <- pars.tv["lambda"]

            Qs_0 <- min(Q[1:10])
            Qq_0 <- 0
            Qst1 <- Qs_0
            Qqt1 <- Qq_0
            cDiffFun <- function(U, Qstar) {
                Qtest <- (1 - alpha_q) * U +
                    v_s_0 * (alpha_q - alpha_s) * U ^ (lambda + 1)
                return(abs(Qstar - Qtest))
            }
            for (t in 1:length(Q)) {
                Qstar <- Q[t] - alpha_s * Qst1 - alpha_q * Qqt1
                ## trivial case: P is zero, therefore U is zero
                if ((P[t] == 0) || (Qstar <= 0)) {
                    U[t] <- 0
                } else {
                    ## estimate U[t] from Q[t], Q[t-1], ...
                    ## we need to know v_s, but this depends on U[t] ^ lambda.
                    ## constraints:
                    v_s_max <- min(1, v_s_0)
                    ## bound based on v_s as function of U (U <= P)
                    v_s_min <- max(0, min(v_s_max, v_s_0 * P[t] ^ lambda))
                                        #v_s_max <- max(v_s_max, v_s_min)
                    beta_max <- v_s_min * (1 - alpha_s) + (1 - v_s_min) * (1 - alpha_q)
                    beta_min <- v_s_max * (1 - alpha_s) + (1 - v_s_max) * (1 - alpha_q)
                    Umin <- Qstar / beta_max
                    Umax <- Qstar / beta_min
                    ## bound based on rainfall
                    Umax <- min(Umax, P[t])
                    ## bound based on lambda = -1
                    Umin <- max(Umin, (Qstar - v_s_0 * (alpha_q - alpha_s)) /
                                 (1 - alpha_q))
                    Umin <- min(Umin, Umax)
                    ## ensure that U values imply possible values of v_s
                    if ((v_s_0 * Umax ^ lambda) >= 1)
                        Umin <- Umax
                    ## check whether we have constrained U to a unique solution
                    if ((Umax - Umin) < 0.0001) {
                        U[t] <- Umax
                    } else {
                        ## approximation based on linear interpolation - any use?
                        #Uest <- Umin + (Umax - Umin) * (1+lambda)
                        result <- optimise(cDiffFun, interval = c(Umin, Umax), Qstar = Qstar)
                        U[t] <- result$minimum
                    }
                }
                ## calculate flow components
                v_s_t <- 0
                if (U[t] > 0)
                    v_s_t <- v_s_0 * (U[t] ^ lambda)
                v_s_t <- max(0, min(1, v_s_t))
                beta_s_t <- v_s_t * (1 - alpha_s)
                beta_q_t <- (1 - v_s_t) * (1 - alpha_q)
                Qst1 <- alpha_s * Qst1 + beta_s_t * U[t]
                Qqt1 <- alpha_q * Qqt1 + beta_q_t * U[t]
            }
        } else {
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
        }
    } else {
        if ("lambda" %in% parSymbols) {
            pars.tv <- tfParsConvert(pars, "tau,v")
            Qs_0 <- min(Q[1:10])
            Qq_0 <- 0
            U <- .C(inverse_filter_lambda,
                    as.double(Q),
                    as.integer(length(Q)),
                    as.double(Qs_0),
                    as.double(Qq_0),
                    as.double(pars.tv["tau_s"]),
                    as.double(pars.tv["tau_q"]),
                    as.double(pars.tv["v_s"]),
                    as.double(pars.tv["lambda"]),
                    as.double(P),
                    U = double(length(Q)),
                    DUP=FALSE, PACKAGE="hydromad")$U
        } else {
            U <- .C(inverse_filter,
                    as.double(Q),
                    as.integer(length(Q)),
                    as.double(a),
                    as.double(b),
                    as.integer(n),
                    as.integer(m),
                    as.double(P),
                    U = double(length(Q)),
                    DUP=FALSE, PACKAGE="hydromad")$U
        }
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
    ## tf.inverse.sim should be in sync with P not Q,
    ## so delay still applies
    if (delay > 0) {
        U <- shiftWindow(U, delay = -delay)
        ## (and fix up time series attributes)
        U <- lag(U, delay)
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
