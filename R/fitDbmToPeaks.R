## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

fitDbmToPeaks <-
    function(MODEL,
             objective = hydromad.getOption("objective"),
             P.quantile = 0.9,
             power = c(0.05, 0.95), qlag = -2:2,
             delay = hydromad.getOption("delay"),
             return_fit = FALSE)
{
    stopifnot(identical(MODEL$sma, "dbm"))
    if (is.na(delay))
        delay <- estimateDelay(MODEL$data[,c("P","Q")])
    P <- observed(MODEL, item = "P", all = TRUE)
    P.thresh <- quantile(coredata(P), P.quantile, na.rm = TRUE)
    dat <- MODEL$data
    dat$P[P < P.thresh] <- NA
    peakMod <- update(MODEL, newdata = dat, power = power, qlag = qlag,
                      routing = "armax", delay = delay, rfit = NULL)
    bestMod <- MODEL
    bestObjVal <- Inf
    lapply(qlag, function(qlagi) {
        modi <- fitByOptim1(update(peakMod, qlag = qlagi), objective = objective)
        obji <- objFunVal(modi, objective = objective)
        if (hydromad.getOption("trace"))
            message("qlag = ", qlagi,
                    "; power = ", signif(coef(modi)[["power"]], 3),
                    "; obj.val = ", signif(obji, 3))
        if (obji < bestObjVal) {
            if (return_fit) {
                bestMod <<- modi
            } else {
                bestMod <<- update(MODEL, power = coef(modi)[["power"]], qlag = qlagi,
                                   and.rescale = FALSE)
            }
            bestObjVal <<- obji
        }
    })
    bestMod
}



## not used:

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

