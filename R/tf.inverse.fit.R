## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

tf.lambda.inverse.fit <-
    function(DATA,
             ..., fit.method = "lambda") #with.lambda = TRUE)
{
    obj <- tf.inverse.fit(DATA, ..., fit.method = fit.method)#with.lambda = with.lambda)
    obj$call <- match.call()
    obj
}

tf.fft.inverse.fit <-
    function(DATA,
             ..., fft.inverse.sim = TRUE)
{
    obj <- tf.inverse.fit(DATA, ..., fft.inverse.sim = fft.inverse.sim)
    obj$call <- match.call()
    obj
}

tf.inverse.fit <-
    function(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             normalise = hydromad.getOption("normalise"),
             fit.method = hydromad.getOption("inverse.fit.method"),
             init.U = TRUE,
             pars = NULL,
             use.Qm = TRUE,
             fft.inverse.sim = FALSE,
             rises.only = FALSE,
             ...,
             max.iterations = hydromad.getOption("inverse.iterations"),
             rel.tolerance = hydromad.getOption("inverse.tolerance"),
             par.epsilon = hydromad.getOption("inverse.epsilon"),
             init.attempt = 0,
             trace = hydromad.getOption("trace"))
{
    DATA <- as.ts(DATA)
    if ("Q" %in% colnames(DATA)) {
        Q <- DATA[,"Q"]
    } else if (NCOL(DATA) == 1) {
        Q <- DATA
    } else {
        stop("need Q in DATA")
    }
    P <- if ("P" %in% colnames(DATA)) DATA[,"P"]

    if (is.na(delay)) {
        delay <- 0
        if ("P" %in% colnames(DATA))
            delay <- estimateDelay(DATA, plot = FALSE)
    }
    with.lambda <- any(grepl("lambda", fit.method)) #(fit.method == "lambda")
    if (with.lambda)
        order <- c(n = 2, m = 1)
    if (!is.null(pars)) {
        pars <- tfParsConvert(pars, "a,b")
        tfParsCheck(pars)
    }
    ## set NAs in flow to zero (insert them again afterwards)
    isna <- is.na(Q)
    Q[isna] <- 0
    sumQ <- sum(Q)

    U <- NULL
    ## option to initialise from U (as rises(Q)), rather than pars.init
    if (init.U) {
        if (!is.null(P)) {
            zeros <- mean(zapsmall(P, 2) == 0, na.rm = TRUE)
        } else {
            zeros <- mean(zapsmall(diff(Q), 2) <= 0, na.rm = TRUE)
        }
        ## backwards difference, i.e. rises are placed at their end time
        rises <- function(x) { x[] <- c(0, pmax(diff(x), 0)); x }
        ## estimate U as the rises in Q, scaled
#        U <- rises(Q)
#        U[U <= quantile(U, zeros, na.rm=TRUE)] <- 0
        ## scale to ensure mass balance with Q
#        U <- U * sum(Q[!isna]) / sum(U[!isna])
        if (!is.null(P)) {
            U <- P
        } else {
            U <- rises(Q)
        }

        ## estimate U as scaled P, scaled in a moving window
        ## (runoff coefficient)
        scale.window <- min(autocorrTime(Q) * 1.5, 10)
        sc <- easysmooth(cbind(Q, U), width = scale.window)
        sc <- sc[,"Q"] / sc[,"U"]
        sc <- na.locf(na.locf(sc, na.rm = FALSE), fromLast = TRUE, na.rm = FALSE)
        sc[!is.finite(sc)] <- 0
        U <- U * sc

    } else {
        ## generate starting parameters
        if (is.null(pars)) {
            pars <- tf.pars.init(DATA, order = order, delay = delay,
                                      with.lambda = with.lambda,
                                      init.attempt = init.attempt)
        }
        ## TODO: hydromad.getOption("catch.errors")
        pcheck <- try(tfParsCheck(pars))
        if (!isTRUE(pcheck))
            return(pcheck)
    }
    oldU <- NULL
    oldpars <- pars
    i <- 1
    repeat {
        if (init.U && (i == 1)) {
            ## we already have an estimate for U
        } else {
            U <- tf.inverse.sim(Q, P = P,
                                     pars = pars, delay = delay,
                                     use.Qm = use.Qm,
                                     use.fft.method = fft.inverse.sim,
                                     rises.only = rises.only)
        }
        ## tf.inverse.sim should be in sync with P not Q,
        ## so delay still applies
        if (with.lambda && !is.null(pars)) {
            mod <- tf.fit(cbind(Q = Q, U = U), method = fit.method,
                               order = order, delay = delay,
                               lambda.init = pars[["lambda"]], ...)
        } else {
            mod <- tf.fit(cbind(Q = Q, U = U), method = fit.method,
                               order = order, delay = delay, normalise = FALSE,
                               ...)
        }
        if (!inherits(mod, "tf"))
            return(mod)

        pars <- coef(mod)

        if (trace) {
            message(paste(" iteration = ", i))
            nicePars <- pars
            if (isStandardModelOrder(pars))
                nicePars <- coef(mod, "tau,v")
            print(nicePars)
        }

        if (!is.null(oldU)) {
            delta <- sum(abs(U - oldU), na.rm = TRUE) / sumQ
            if (delta < rel.tolerance) break
        }

        if (FALSE && !is.null(oldpars)) {
            ## select and order by name for direct comparison
            pars <- pars[names(oldpars)]
            names(pars) <- names(oldpars)
            ## omitted parameter values are assumed to be implicitly 0
            pars[is.na(pars)] <- 0
            par.deltas <- abs((pars - oldpars) / oldpars)
            if (trace) print(par.deltas)
            if (max(par.deltas) < par.epsilon) break
        }

        if (i >= max.iterations) {
            warning("tf.inverse.fit reached maximum number of iterations")
            break
        }
        oldU <- U
        oldpars <- pars
        i <- i + 1
    }

    if (normalise)
        mod <- normalise.tf(mod)
    mod$call <- match.call()
    mod
}


tf.pars.init <-
    function(DATA, order = NA, delay = NA,
             normalise = hydromad.getOption("normalise"),
             with.lambda = FALSE, init.attempt = 0)
{
    if (with.lambda)
        order <- c(n = 2, m = 1)
    stopifnot(length(order) == 2)
    if (is.null(names(order)))
        names(order) <- c("n", "m")
    n <- order["n"]
    m <- order["m"]
    #alphas <- alpha.init(n = order["n"])

    stdOrder <- (identical(order, c(n=0,m=0)) ||
                 identical(order, c(n=1,m=0)) ||
                 identical(order, c(n=1,m=1)) ||
                 identical(order, c(n=2,m=1)) ||
                 identical(order, c(n=2,m=2)) ||
                 (with.lambda))

    if (!stdOrder) {
        Q <- coredata(DATA[,"Q"])
        theta <- coef(arima(Q, c(order[1], 0, order[2]),
                            include.mean = FALSE))
        ## arima has an implicit value of 1 for b_0
        theta <- c(theta[seq_len(n)],
                   1, theta[n + seq_len(m)])
        if (normalise)
            theta <- normalise.tf.coef(theta)
        names(theta) <-
            c(if (n > 0)
              paste("a", seq_len(n), sep="_"),
              paste("b", 0:m, sep = "_"))
        return(theta)
    }

    set.seed(init.attempt)
    pctile.s <- runif(1, min=0.1, max = 1) # 0.5
    pctile.q <- runif(1, min=0, max = 0.1) # 0.02
    v_s <- runif(1, min = 0, max = 1) # 0.5
    v_3 <- runif(1, min = 0, max = (1 - v_s))
    v_q <- (1 - v_s)
    if (order["m"] >= 2)
        v_q <- (1 - v_s - v_3)
    lambda <- - runif(1, min = 0, max = 0.9) # -0.01

    Q <- coredata(DATA[,"Q"])
    qrecc <- c(NA, Q) / c(Q, NA)
    qrecc <- qrecc[is.finite(qrecc)]
    qrecc <- qrecc[qrecc < 1]
    alpha_s <- quantile(qrecc, pctile.s, names = FALSE, na.rm = TRUE)
    alpha_q <- quantile(qrecc, pctile.q, names = FALSE, na.rm = TRUE)
    tau_s <- -1/log(alpha_s)
    tau_q <- -1/log(alpha_q)

    pars <-
     c(if (order["n"] >= 1) c(tau_s = tau_s),
      if (order["n"] >= 2) c(tau_q = tau_q),
      if (order["m"] >= 1) c(v_s = v_s),
      if (order["m"] >= 1) if (!with.lambda) c(v_q = v_q),
      if (order["m"] >= 2) c(v_3 = v_3),
      if (with.lambda) c(lambda = lambda))

    tfParsConvert(pars, "a,b")
}


