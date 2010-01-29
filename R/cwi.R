## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## Catchment Wetness Index (CWI) soil moisture accounting
## aka "classic"
## with threshold parameter extensions of Ye et. al. (1997)
## tw = drying rate at reference temperature
## f = temperature dependence of drying rate
## c = mass balance term
## l = moisture threshold for producing flow
## p = power on soil moisture
## t_ref = reference temperature (traditionally 20 deg C)
cwi.sim <-
    function(DATA,
             tw, f = 0, c,
             l = 0, p = 1,
             t_ref = hydromad.getOption("cwi")$t_ref,
             s_0 = 0,
             return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    if (NCOL(DATA) > 1) stopifnot("P" %in% colnames(DATA))
    if (f > 0) stopifnot("E" %in% colnames(DATA))
    P <- if (NCOL(DATA) > 1) DATA[,"P"] else DATA
    ## special value c = NA used for initial run for scaling
    if (is.na(c))
        c <- 1
    ## check values
    stopifnot(tw >= 0)
    stopifnot(f >= 0)
    stopifnot(c >= 0)
    stopifnot(l >= 0)
    stopifnot(p > 0)
    stopifnot(length(t_ref) == 1)
    ## compute soil moisture index s
    if (f == 0) {
        ## this is an invariant drying rate model
        w <- 1 - 1 / max(tw, 1)
        w <- rep(w, length(P))
        ## use filter_tv rather than filter because it maintains
        ## state over NA gaps (filter resets to 0)
        #s <- filter(P, w, method="recursive", init = s_0)
    } else {
        tw_k <- tw * exp(- 0.062 * f * (DATA[,"E"] - t_ref))
        tw_k <- pmax(tw_k, 1)
        w <- 1 - 1 / tw_k
    }
    s <- filter_tv(P, w, init = s_0)
    ## compute effective rainfall U
    U <- (pmax(s - l, 0) ^ p) * P
    ## TODO: test mass balance scaling with p > 1
    U <- (c ^ p) * U
    if (return_state) return(ts.union(U=U, s=s, w=w))
    return(U)
}


absorbScale.cwi <- function(object, gain)
{
    #if (!isTRUE(object$estimateScale)) {
    coeff <- coef(object)
    c <- coeff[["c"]]
    ## we only want to do this when c is NA (special value)
    if (is.null(c) || !is.na(c))
        return(NULL)
    c <- 1
    p <- coef(object)[["p"]]
    c <- c * (gain ^ p)
    c <- max(c, 0)
    object$parlist[["c"]] <- c
    object$U <- (c ^ p) * object$U
    object
}
