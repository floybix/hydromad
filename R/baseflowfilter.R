## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

baseflowfilter <- function(x, ...)
    UseMethod("baseflowfilter")

baseflowfilter.hydromad <- function(x, ...)
{
    ## TODO
    armax.inverse.sim()
}

baseflowfilter.default <- function(x, ..., method = "rollmin", width)
{
    rollapply(rollapply(x, min, width = width, ...), mean, width = width, ...)
}
