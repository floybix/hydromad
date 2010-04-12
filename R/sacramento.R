## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## Sacramento Soil Moisture Accounting model.
## Developed by the US National Weather Service.
sacramento.sim <-
    function(DATA,
             uztwm, uzfwm, uzk, pctim, adimp, zperc, rexp,
             lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree,
             etmult = 1, dt = 1,
             return_state = FALSE)
{
    stopifnot(c("P","E") %in% colnames(DATA))
    ## check values
    stopifnot(uztwm >= 0)
    stopifnot(uzfwm >= 0)
    stopifnot(uzk >= 0)
    stopifnot(0 <= pctim && pctim <= 1)
    stopifnot(adimp >= 0)
    stopifnot(zperc >= 0)
    stopifnot(lztwm >= 0)
    stopifnot(lzfsm >= 0)
    stopifnot(lzfpm >= 0)
    stopifnot(lzsk >= 0)
    stopifnot(lzpk >= 0)
    stopifnot(pfree >= 0)
    stopifnot(etmult >= 0)
    stopifnot(dt >= 0)

    xpar <-
        c(uztwm, uzfwm, uzk, pctim, adimp, zperc, rexp,
          lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree)

    P <- DATA[,"P"]
    E <- DATA[,"E"]
    ## skip over missing values
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    ## NOTE: there is not ability to return state variables yet
    ## also no pure.R.code implementation
    U <- .C(sma_sac,
            as.double(P),
            as.double(E),
            as.integer(NROW(DATA)),
            as.double(xpar),
            as.double(etmult),
            as.double(dt),
            U = double(NROW(DATA)),
            NAOK = FALSE, DUP = FALSE, PACKAGE="hydromad")$U
    ## make it a time series object again
    attributes(U) <- attributes(P)
    ## re-insert missing values
    U[bad] <- NA
    return(U)
}

sacramento.ranges <- function()
    list(uztwm = c(1, 150),
         uzfwm = c(1, 150),
         uzk = c(0.1, 0.5),
         pctim = c(0.000001, 0.1),
         adimp = c(0, 0.4),
         zperc = c(1, 250),
         rexp = c(0, 5),
         lztwm = c(1, 500),
         lzfsm = c(1, 1000),
         lzfpm = c(1, 1000),
         lzsk = c(0.01, 0.25),
         lzpk = c(0.0001, 0.25),
         pfree = c(0, 0.6))
