## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##

## Catchment Moisture Deficit (CMD) soil moisture accounting
## d = CMD flow threshold
## f = CMD stress threshold
## e = temperature to PET conversion factor
cmd.sim <-
    function(DATA,
             d, f, e, M_0 = d/2,
             return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    stopifnot(c("P","E") %in% colnames(DATA))
    ## check values
    stopifnot(d >= 0)
    stopifnot(f >= 0)
    stopifnot(e >= 0)
    ## default initial state
    if (is.na(M_0)) M_0 <- d / 2

    ## f is expressed as a proportion of d
    g <- f * d

    P <- DATA[,"P"]
    E <- DATA[,"E"]
    ## skip over missing values (maintaining the state M)
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    ## TODO: return state from C code
    COMPILED <- (ihacres.getOption("pure.R.code") == FALSE)
    if (COMPILED && !return_state) {
        U <- .C(sma_cmd,
                as.double(P),
                as.double(E),
                as.integer(NROW(DATA)),
                as.double(d),
                as.double(g),
                as.double(e),
                as.double(M_0),
                U = double(NROW(DATA)),
                NAOK=FALSE, DUP=FALSE, PACKAGE="ihacreslab")$U
        ## make it a time series object again
        mostattributes(U) <- attributes(DATA)
        class(U) <- "ts"
    } else {
        ## implementation in R for cross-checking (slow)
        U <- M <- ET <- P
        M_prev <- M_0
        for (t in seq(1, length(P))) {
            ## rainfall reduces CMD (Mf)
            if (P[t] > 0) {
                if (M_prev < d)
                    Mf <- M_prev * exp(-P[t] / d)
                else if (M_prev < d + P[t])
                    Mf <- d * exp((-P[t] + M_prev - d) / d)
                else
                    Mf <- M_prev - P[t]
            } else {
                Mf <- M_prev
            }
            ## drainage (rainfall not accounted for in -dM)
            U[t] <- max(0, P[t] - M_prev + Mf)
            ## evapo-transpiration
            ET[t] <- e * E[t] * min(1, exp(2 * (1 - Mf / g)))
            ET[t] <- max(0, ET[t])
            ## mass balance
            M[t] <- M_prev - P[t] + U[t] + ET[t]
            M_prev <- M[t] <- max(0, M[t])
            ##M[t] <- M_prev <- max(0, Mf + ET[t])
        }
    }
    ## re-insert missing values
    U[bad] <- NA
    if (return_state) return(ts.union(U=U, CMD=M, ET=ET))
    return(U)
}
