## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## the linear flow routing module
## (calculates streamflow from effective rainfall)
## tf.sim(U, c(v_s))                    ## scaling of input, no dynamics
## tf.sim(U, c(tau_s))                  ## single exponential store
## tf.sim(U, c(tau_s, v_s))             ## exponential store and instantaneous store in parallel (v_q = 1 - v_s)
## tf.sim(U, c(tau_s, tau_q, v_s))      ## two exponential stores in parallel
## tf.sim(U, c(tau_s, tau_q, v_s, v_3)) ## two exponential stores and instantaneous store in parallel (v_q = 1 - v_s - v_3)
## tf.sim(U, c(tau_s, tau_q, tau_3, v_s, v_3)) ## three exponential stores in parallel (v_q = 1 - v_s - v_3)
## tf.sim(U, c(tau_s, tau_q, v_s, lambda)) ## two exponential stores in parallel, with v_s = v_s * (U ^ lambda)
tf.sim <-
    function(U,
             pars = c(tau_s=0, tau_q=0, v_s=1, v_3=0),
             delay = 0, Xs_0 = 0, Xq_0 = 0, X3_0 = 0,
             return_components = FALSE,
             na.action = na.pass,
             epsilon = hydromad.getOption("sim.epsilon"))
{
    ## note U is allowed to be multi-variate, i.e. multiple columns
    if (!is.ts(U)) U <- as.ts(U)
    U <- na.action(U)
    ## apply 'U' delay in reverse to 'X' (i.e. lag by 'delay' steps)
    ## so can take delay as 0 for simulation purposes
    if (delay != 0)
        U <- lag(U, -delay)
    ## null case of no parameters
    if (length(pars) == 0) {
        X <- U
        X <- shiftWindow(X, delay)
        X[X < epsilon] <- 0
        return(X)
    }
    pars <- unlist(pars)
    ## check parameters in ARMAX form, even though not used directly here
    ## could just be a warning?
    tfParsCheck(pars)
    ## extract first parts of parameter names (before underscore)
    parSymbols <- gsub("_.*", "", names(pars))

    ## check if given parameters are like a_1, a_2, b_0, b_1...
    if (any(parSymbols %in% c("a", "b"))) {
        ## ARMAX model formulation
        a <- pars[parSymbols == "a"]
        b <- pars[parSymbols == "b"]
        stopifnot(length(b) > 0)
        n <- length(a)
        m <- length(b) - 1

        if ((return_components == FALSE) &&
            (("lambda" %in% parSymbols) == FALSE) &&
            (("loss" %in% parSymbols) == FALSE))
        {
            ## run ARMAX model directly
            a_init <- rep(Xs_0, n)
            X <- filter(U, b, sides = 1)
            if (NCOL(U) > 1) {
                X[seq_len(m),] <- U[seq_len(m),] * b[1]
            } else {
                X[seq_len(m)] <- U[seq_len(m)] * b[1]
            }
            if (length(a) > 0)
                X <- filter(X, a, method = "recursive", init = a_init)
            ## zap simulated values smaller than epsilon
            X[X < epsilon] <- 0
            X <- shiftWindow(X, delay)
            return(X)
        }

        ## else... convert to exponential components formulation
    }

    ## otherwise: exponential components formulation
    pars <- tfParsConvert(pars, "tau,v")
    parSymbols <- gsub("_.*", "", names(pars))
    ## check if given parameters are like tau_s, tau_q, v_s, v_q...
    if (!any(parSymbols %in% c("tau", "v")))
        stop("items in 'pars' should have names starting a_ / b_ / tau_ / v_")
    ## fill in default values
    tau_s <- 0
    tau_q <- 0
    tau_3 <- 0
    v_s <- 1
    v_3 <- 0
    series <- 0
    lambda <- 0
    loss <- 0
    with(as.list(pars), {
        if (("v_q" %in% names(pars)) == FALSE)
            v_q <- max(1 - v_s - v_3, 0)
        stopifnot(all(c(tau_s, tau_q, Xs_0, Xq_0) >= 0))
        ## convert from 'tau' and 'v' to 'alpha' and 'beta'
        alpha_s <- exp(-1 / tau_s)
        alpha_q <- exp(-1 / tau_q)
        alpha_3 <- exp(-1 / tau_3)
        ## lambda parameter defines dependence of v_s on U
        if (lambda != 0) {
            v_s <- v_s * (U ^ lambda)
            v_s <- pmax(0, pmin(1, v_s)) ## ensure (0 <= v_s <= 1)
            v_q <- pmax(0, 1 - v_s - v_3)
        }
        ## note: here v_s, v_q can be either scalars or vectors!
        beta_s <- v_s * (1 - alpha_s)
        beta_q <- v_q * (1 - alpha_q)
        beta_3 <- v_3 * (1 - alpha_3)
        ## apply exponential decay filter to each component
        ## note filter_loss is equivalent to filter if loss = 0
        ## convert loss from G[k] model into a Q[k] formulation
        lossVal <- (1 - alpha_s) * loss
        Xs <- filter_loss(beta_s * U, alpha_s, loss = lossVal, init = Xs_0)
        Xq <- filter(beta_q * U, alpha_q, method = "recursive", init = Xq_0)
        X3 <- if (v_3) filter(beta_3 * U, alpha_3, method = "recursive", init = X3_0)

        ## align results to original input
        Xs <- shiftWindow(Xs, delay)
        Xq <- shiftWindow(Xq, delay)
        X3 <- if (v_3) shiftWindow(X3, delay)

        ## zap simulated values smaller than epsilon
        Xs[Xs < epsilon] <- 0
        Xq[Xq < epsilon] <- 0
        if (v_3) X3[X3 < epsilon] <- 0

        ## can only simulate stores in parallel here (TODO)
        if (series > 0)
            warning("simulation not correct: series=", series,
                    " but assumed stores in parallel")

        if (return_components) {
            if (v_3) return(ts.union(Xs = Xs, Xq = Xq))
            else return(ts.union(Xs = Xs, Xq = Xq, X3 = X3))
        } else {
            if (v_3) return(Xs + Xq + X3)
            else return(Xs + Xq)
        }
    })
}

