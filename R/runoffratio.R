## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

runoffratio.sim <-
    function(DATA, width = NA, c = 2, sides = 2, qlag = 0, return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    stopifnot(c("P", "Q") %in% colnames(DATA))
    P <- DATA[,"P"]
    Q <- DATA[,"Q"]
    ## check values
    stopifnot(c >= 1)
    ## synchronise Q lagged by qlag
    Qd <- shiftWindow(Q, round(qlag), and.lag = TRUE)
    ## compute effective rainfall U
    ## estimate U as scaled P, scaled in a moving window
    ## (runoff coefficient)
    if (is.na(width)) {
        width <- autocorrTime(Q)
        width <- min(max(width, 8), length(Q) %/% 5)
    }
    sc <- simpleSmoothTs(cbind(Q, P), width = width, c = c, sides = 2)
    sc <- sc[,"Q"] / sc[,"P"]
    sc <- na.locf(na.locf(sc, na.rm = FALSE), fromLast = TRUE, na.rm = FALSE)
    sc[!is.finite(sc)] <- 0
    U <- P * sc
    U
}
