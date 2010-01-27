## ihacres: IHACRES rainfall runoff model
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
## IHACRES model by Jakeman, Croke and others
##

flow.inputs.sim <-
    function(n = 600,
             U = ts(pmax(0, rgamma(n, shape=0.1, scale=20) - 5)),
             pars = c(tau_s=64, tau_q=1, v_s=0.5),
             threshq = 0.9, inter = 2,
             ...)
{
    ## normalise transfer function
    Q <- tf.sim(U, pars=pars)
    Qhat <- flow.error.sim(Q) ## TODO: what type of errors?
    Uhat <- flow.inputs(Qhat, method="tf", ...)
    Uhatr <- flow.inputs(Qhat, method="rises", ...)
    thresh <- quantile(Uhat, threshq)
    events <- event.clusters(cbind(U,Uhat), thresh=thresh, inter=inter)
    Usums <- eventapply(cbind(U,Uhat), thresh=thresh, inter=inter, FUN=sum)
    Usumsr <- eventapply(cbind(U,Uhatr), thresh=thresh, inter=inter, FUN=sum)
    list(U=U, Uhat=Uhat, Uhatr=Uhatr, Q=Q, Qhat=Qhat, Usums=Usums, Usumsr=Usumsr, events=events)
}

rain.error.sim <-
    function(P,
             err.frac = 0.3,
             resamp.frac = 0.1,
             lags = c(rep(1, n*0.3), rep(0, n*0.7)),
             ...)
{
    n <- length(P)
    ## simulate error in areal rain by scaling randomly...
    P.hat <- P * runif(n, min=1-err.frac, max=1+err.frac)
    ## random re-sampling 10 pct of days...
    nr <- n * resamp.frac
    P.hat[sample(seq(along=P), nr)] <- sample(P, nr)
    ## ...and scaling a random subset (of half the data) by half
                                        #ii <- sample(1:length(P), length(P)/2)
                                        #P.hat[ii] <- P.hat[ii] / 2
    ## timing error: lag first third of data
    ## TODO lags
    P.hat[1:100] <- c(NA, P.hat[1:99])
}

flow.error.sim <-
    function(x,
             add.err.frac = 0.01,
             mult.err.frac = 0.1,
             ma.window = 20,
             sd.mult.err = 0.1,
             scale.trim = 0.1,
             digits.zap = 3)
{

    ## TODO: rating curve error

    xscale <- mean(x, trim = scale.trim, na.rm = TRUE)
    sd.mult.err <- mult.err.frac * xscale
    sd.add.err <- add.err.frac * xscale
    ## apply multiplicative error to observed flow
    ## uniform random, averaged over a moving window (i.e. auto-correlated)
    err <- runif(length(x))
    err <- filter(err, rep(1, ma.window), circular=TRUE)
    err <- ts(pmax(0, 1 + scale(err) * sd.mult.err))
    ## uncorrelated additive error
    adderr <- rnorm(length(x), sd = sd.add.err)
    xhat <- pmax(x * err + adderr, 0)
    ## rounding error
    xhat <- zapsmall(xhat, digits = digits.zap)
    xhat
}
