## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## a = time to drop to half
powuh.sim <-
    function(U, delay = 0,
             a, b = 1, c = 1,
             init = 0, na.action = na.pass,
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

    ## TODO

    t <- 0:100
    uh <- (1 / (1 + (t/a)^(b/c))) ^ c

    X <- filter(U, uh, sides = 1)
    
    ## align results to original input
    X <- shiftWindow(X, delay)

    ## zap simulated values smaller than epsilon
    X[X < epsilon] <- 0
    X
}


ssg.powuh <- function(theta)
{
    if (length(theta) == 0)
        return(1)
    1
    #theta <- tfParsConvert(theta, "a,b")
    #ssg.tf.coef(theta)
}

normalise.powuh <- function(theta)
{
    theta
}
