## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

scalar.sim <-
    function(DATA, frac, return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    if (NCOL(DATA) > 1) stopifnot("P" %in% colnames(DATA))
    P <- if (NCOL(DATA) > 1) DATA[,"P"] else DATA
    ## special value frac = NA used for initial run for scaling
    if (is.na(frac))
        frac <- 1
    ## check values
    stopifnot(frac >= 0)
    ## compute effective rainfall U
    frac * P
}

scalar.ranges <- function()
    list(frac = NA)

absorbScale.scalar <- function(object, gain)
{
    coeff <- coef(object)
    frac <- coeff[["frac"]]
    ## we only want to do this when frac is NA (special value)
    if (is.null(frac) || !is.na(frac))
        return(NULL)
    frac <- 1
    frac <- frac * gain
    object$parlist[["frac"]] <- frac
    object$U <- frac * object$U
    object
}
