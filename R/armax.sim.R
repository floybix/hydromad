## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

armax.sim <-
    function(U, a_1 = 0, a_2 = 0, a_3 = 0, 
             b_0 = 1, b_1 = 0, b_2 = 0, b_3 = 0,
             pars = NULL,
             delay = 0, init = 0, na.action = na.pass,
             epsilon = hydromad.getOption("sim.epsilon"))
{
    ## parameter vectors
    a <- c(a_1, a_2, a_3)
    b <- c(b_0, b_1, b_2, b_3)
    a <- stripzeros(a)
    b <- stripzeros(b, up.to = 1)
    ## allow parameters to be passed in a single named vector
    if (length(pars) > 0) {
        pars <- tfParsConvert(pars, "a,b")
        a <- pars[grep("^a", names(pars))]
        b <- pars[grep("^b", names(pars))]
        stopifnot(length(b) > 0)
    }
    ## model order
    n <- length(a)
    m <- length(b) - 1
    ## check values
    stopifnot(all(is.finite(c(a, b))))
    ## stability check (based on stats::arima)
    if (length(a) > 0)
        if (!all(Mod(polyroot(c(1, -a))) > .95))
            stop("AR component not in the region of stationarity")
    ## note U is allowed to be multi-variate, i.e. multiple columns
    if (!is.ts(U)) U <- as.ts(U)
    U <- na.action(U)
    ## apply 'U' delay in reverse to 'X' (i.e. lag by 'delay' steps)
    ## so can take delay as 0 for simulation purposes
    if (delay != 0)
        U <- lag(U, -delay)
    ## run ARMAX model
    init <- rep(init, length = n)
    X <- filter(U, b, sides = 1)
    if (length(a) > 0) {
        X <- shiftWindow(X, -m)
        X <- filter(X, a, method = "recursive", init = init)
        X <- shiftWindow(X, m)
    }
    ## zap simulated values smaller than epsilon
    X[abs(X) < epsilon] <- 0
    X <- shiftWindow(X, delay)
    X
}


ssg.armax <- function(theta)
    ssg.tf.coef(theta)

normalise.armax <- function(theta)
    normalise.tf.coef(theta)

