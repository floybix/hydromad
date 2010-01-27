## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##

## from Farmer et al 2003, WRR
## general mass balance structure:
## dS/dt = p - q(S) - e(S, Ep)
## slightly modified forms as given in Bai et al 2009, EMS


bucket.sim <-
    function(DATA,
             Sb, fc, a.ei = 0, M = 1, a.ss = 0, S_0 = 0,
             return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    stopifnot(c("P","E") %in% colnames(DATA))
    ## check values
    stopifnot(Sb >= 0)
    stopifnot(S_0 >= 0)
    stopifnot(0 <= fc && fc <= 1)
    stopifnot(0 <= a.ei && a.ei <= 1)
    stopifnot(0 <= M && M <= 1)
    stopifnot(0 <= a.ss && a.ss <= 1)

    ## fc is expressed as a proportion of Sb
    Sfc <- fc * Sb

    P <- DATA[,"P"]
    E <- DATA[,"E"]
    ## skip over missing values (maintaining the state S)
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    ## TODO: return state from C code
    COMPILED <- (ihacres.getOption("pure.R.code") == FALSE)
    if (FALSE && COMPILED && !return_state) {
        U <- .C(NA, #sma_bucket, TODO
                as.double(P),
                as.double(E),
                as.integer(NROW(DATA)),
                as.double(Sb),
                as.double(fc),
                as.double(a.ei),
                as.double(M),
                as.double(a.ss),
                as.double(S_0),
                U = double(NROW(DATA)),
                NAOK=FALSE, DUP=FALSE, PACKAGE="ihacreslab")$U
        ## make it a time series object again
        mostattributes(U) <- attributes(DATA)
        class(U) <- "ts"
    } else {
        ## implementation in R for cross-checking (slow)
        U <- S <- ET <- P
        S_prev <- S_0
        for (t in seq(1, length(P))) {
            ## evapo-transpiration
            Eintc <- a.ei * P[t]
            S[t] <- min(Sb, S_prev + P[t] - Eintc)
            Etrans <- M * min(1, S[t] / Sfc) * E[t]
            Ebare <- (1 - M) * (S[t] / Sb) * E[t]
            ET[t] <- Eintc + Etrans + Ebare
            ## mass balance
            S[t] <- S_prev + P[t] - ET[t]
            ## drainage (saturation excess)
            Use <- max(0, S[t] - Sb)
            S[t] <- S[t] - Use
            ## drainage (sub-surface)
            Uss <- a.ss * max(0, S[t] - Sfc)
            S[t] <- S[t] - Uss
            U[t] <- Use + Uss
            S_prev <- S[t]
        }
    }
    ## re-insert missing values
    U[bad] <- NA
    if (return_state) return(ts.union(U=U, S=S, ET=ET))
    return(U)
}
