## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## a = time to drop to half
powuh.sim <-
    function(U, delay = 0,
             a = 1, b = 1, c = 1, init = 0, uhsteps = 100,
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

    t <- seq(0, uhsteps)
    uh <- (1 / (1 + (t/a)^(b/c))) ^ c
    ## normalise:
    uh <- uh / sum(uh)
    ## initialisation
    Upad <- window(U, start = start(U)[1] - uhsteps, extend = TRUE)
    Upad[1:uhsteps] <- init
    X <- filter(U, uh, sides = 1)
    X <- window(X, start = start(U))
    
    ## align results to original input
    X <- shiftWindow(X, delay)

    ## zap simulated values smaller than epsilon
    X[X < epsilon] <- 0
    X
}


ssg.powuh <- function(theta)
{
    1
}

normalise.powuh <- function(theta)
{
    theta
}
