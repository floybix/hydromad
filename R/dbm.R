## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

fitDbmByGam <-
    function(MODEL, ...)
{
    start_time <- proc.time()
    parlist <- as.list(coef(MODEL, warn = FALSE))
    qlag <- parlist$qlag
    ## TODO: test each qlag
    qlag <- 0
    Q <- MODEL$data[,"Q"]
    P <- MODEL$data[,"P"]
    ## TODO: should use a time-varying parameter model?
    pqTf <- armax.sriv.fit(cbind(U = P, Q = Q), order = c(1,0))
    tfc <- coef(pqTf)
    ar1 <- tfc[["a_1"]]
    bb <- pmax(Q - ar1 * lag(Q,-1), 0) / lag(P, qlag)
    ## TODO: use gam to find form of bb ~ Q
    library(mgcv)
    foo <- gam(bb ~ s(Q[-1]), subset = is.finite(bb))
    plot(foo)
    ## TODO: parameterise as power law
    bestModel$funevals <- 0
    bestModel$timing <- proc.time() - start_time
    return(bestModel)

}

dbm.sim <-
    function(DATA, power, c = 1, qlag = 0, return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    stopifnot(c("P", "Q") %in% colnames(DATA))
    P <- DATA[,"P"]
    Q <- DATA[,"Q"]
    ## special value c = NA used for initial run for scaling
    if (is.na(c))
        c <- 1
    ## check values
    stopifnot(c >= 0)
    ## synchronise Q lagged by qlag
    Qd <- shiftWindow(Q, round(qlag), and.lag = TRUE)
    ## compute effective rainfall U
    c * P * Qd ^ power
}

dbm.ranges <- function()
    list(power = c(0.01, 0.9),
         qlag = c(-1, 2),
         c = NA)

absorbScale.hydromad.dbm <- function(object, gain, ...)
{
    ## uses 'c' parameter for scaling:
    absorbScale.hydromad.scalar(object, gain, parname = "c")
}
