
HydroTestData <- local({
    len <- 730
    ## regular rainfall impulse
    P <- ts(rep(0, length = len))
    P[seq(10, len-1, by = 20)] <- 5
    ## make one rain event much larger
    P[quantile(which(P > 0), 0.67, type = 1)] <- 20
    ## sine wave for temperature
    E <- ts(15 + 15 * sin(seq(0, 4*pi, length = len)))
    ## flow based on square of rainfall and inverse to temperature
    Q <- filter(0.01 * P^2 * (1 - E / max(E)),
                filter = c(1.4, -0.45), method = "recursive")
    as.zooreg(ts.intersect(P = P, E = E, Q = Q))
})
