## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

scalar.sim <-
    function(DATA, c, return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    if (NCOL(DATA) > 1) stopifnot("P" %in% colnames(DATA))
    P <- if (NCOL(DATA) > 1) DATA[,"P"] else DATA
    ## special value c = NA used for initial run for scaling
    if (is.na(c))
        c <- 1
    ## check values
    stopifnot(c >= 0)
    ## compute effective rainfall U
    c * P
}

scalar.ranges <- function()
    list(c = NA)

absorbScale.hydromad.scalar <- function(object, gain, parname = "c", ...)
{
    if (gain <= 0)
        return(NULL)
    coeff <- coef(object, which = "sma")
    if (parname %in% names(coeff) == FALSE)
        return(NULL)
    c <- coeff[[parname]]
    ## we only want to do this when c is NA (special value)
    if (is.null(c) || !is.na(c))
        return(NULL)
    c <- 1
    c <- c * gain
    object$parlist[[parname]] <- c
    object$U <- c * object$U
    object
}
