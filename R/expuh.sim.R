## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

expuh.sim <-
    function(U, delay = 0,
             tau_s = 0, tau_q = 0, tau_3 = 0,
             v_s = 1, v_q = NA, v_3 = 0, v_4 = 0,
             series = 0,
             loss = 0,
             Xs_0 = 0, Xq_0 = 0, X3_0 = 0,
             return_components = FALSE,
             na.action = na.pass,
             epsilon = hydromad.getOption("sim.epsilon"))
{
    delay <- round(delay)
    ## note U is allowed to be multi-variate, i.e. multiple columns
    if (!is.ts(U)) U <- as.ts(U)
    U <- na.action(U)
    ## apply 'U' delay in reverse to 'X' (i.e. lag by 'delay' steps)
    ## so can take delay as 0 for simulation purposes
    if (delay != 0)
        U <- lag(U, -delay)

    ## TODO: this should depend on 'series'
    if (is.na(v_q))
        v_q <- max(0, min(1, 1 - v_s - v_3))

    stopifnot(all(c(tau_s, tau_q, tau_3) >= 0))
    ## convert from 'tau' and 'v' to 'alpha' and 'beta'
    alpha_s <- exp(-1 / tau_s)
    alpha_q <- exp(-1 / tau_q)
    alpha_3 <- exp(-1 / tau_3)
    beta_s <- v_s * (1 - alpha_s)
    beta_q <- v_q * (1 - alpha_q)
    beta_3 <- v_3 * (1 - alpha_3)
    if ((series == 0) || return_components) {
        ## components in parallel.
        ## apply exponential decay filter to each component
        ## note filter_loss is equivalent to filter if loss = 0
        ## convert loss from G[k] model into a Q[k] formulation
        lossVal <- (1 - alpha_s) * loss
        Xs <- filter_loss(beta_s * U, alpha_s, loss = lossVal, init = Xs_0)
        Xq <- filter(beta_q * U, alpha_q, method = "recursive", init = Xq_0)
        X3 <- if (v_3) filter(beta_3 * U, alpha_3, method = "recursive", init = X3_0)
    } else {
        if (v_3 == 0) {
            ## second-order model
            Xs <- filter(beta_s * U, alpha_s, method = "recursive", init = Xs_0)
            Xq <- filter(Xs * beta_q * U, alpha_q, method = "recursive", init = Xq_0)
            Xs[] <- X3 <- 0
        } else {
            ## third-order model
            Xs <- filter(beta_s * U, alpha_s, method = "recursive", init = Xs_0)
            if (series == 1) {
                ## two components in series and one in parallel
                ## (s & 3 are in series; q in parallel)
                X3 <- filter(beta_3 * U, alpha_3, method = "recursive", init = X3_0)
                Xs <- filter(X3 * beta_s * U, alpha_s, method = "recursive", init = Xs_0)
                X3[] <- 0
                Xq <- filter(beta_q * U, alpha_q, method = "recursive", init = Xq_0)
            } else if (series == 2) {
                ## one component in series with two in parallel
                ## (3 in series; s & q in parallel)
                X3 <- filter(beta_3 * U, alpha_3, method = "recursive", init = X3_0)
                Xs <- filter(X3 * beta_s * U, alpha_s, method = "recursive", init = Xs_0)
                Xq <- filter(X3 * beta_q * U, alpha_q, method = "recursive", init = Xq_0)
                X3[] <- 0
            } else if (series == 3) {
                ## three components in series
                X3 <- filter(beta_3 * U, alpha_3, method = "recursive", init = X3_0)
                Xs <- filter(X3 * beta_s * U, alpha_s, method = "recursive", init = Xs_0)
                X3[] <- 0
                Xq <- filter(Xs * beta_q * U, alpha_q, method = "recursive", init = Xq_0)
                Xs[] <- 0
            } else {
                stop("unrecognised values of 'series': ", series)
            }
        }
    }

    ## align results to original input
    Xs <- shiftWindow(Xs, delay)
    Xq <- shiftWindow(Xq, delay)
    X3 <- if (v_3) shiftWindow(X3, delay)

    ## zap simulated values smaller than epsilon
    Xs[Xs < epsilon] <- 0
    Xq[Xq < epsilon] <- 0
    if (v_3) X3[X3 < epsilon] <- 0

    if (return_components) {
        if (v_3) return(ts.union(Xs = Xs, Xq = Xq, X3 = X3))
        else return(ts.union(Xs = Xs, Xq = Xq))
    } else {
        if (v_3) return(Xs + Xq + X3)
        else return(Xs + Xq)
    }
}


ssg.expuh <- function(theta)
{
    if (length(theta) == 0)
        return(1)
    theta <- tfParsConvert(theta, "a,b")
    ssg.tf.coef(theta)
}

normalise.expuh <- function(theta)
{
#    v_s <- 1
#    v_3 <- 0
#    v_q <- NA
#    if ("v_s" %in% names(theta))
#        v_s <- theta[["v_s"]]
#    if ("v_3" %in% names(theta))
#        v_3 <- theta[["v_3"]]
#    if ("v_q" %in% names(theta))
#        v_q <- theta[["v_q"]]
#    if (is.na(v_q)) {
#        ## unit volume is enforced anyway
#        return(theta)
#    }
#    vv <- v_s + v_q + v_3

    theta <- tfParsConvert(theta, "a,b")
    tmp <- normalise.tf.coef(theta)
    tfParsConvert(tmp, "tau,v")
}
