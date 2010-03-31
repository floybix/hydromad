## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


predict.hydromad <-
    function(object, newdata = NULL,
             which = c("both", "sma", "routing"),
             ...)
{
    which <- match.arg(which)
    if (is.null(newdata)) {
        if (which == "routing") {
            newdata <- object$U
        } else {
            newdata <- object$data
        }
    } else {
        newdata <- as.ts(newdata)
    }
    ## check that parameters are fully specified
    if (!isFullySpecified(object, which = which))
        stop("model parameters are not fully specified")
    ## construct call to simulation function
    simpars <- coef(object, which = which, etc = TRUE)
    ucall <- quote(hydromad.sim(newdata))
    ucall <- as.call(c(as.list(ucall),
                       as.list(simpars),
                       ...))
    if (which %in% c("sma", "both"))
        ucall["sma"] <- list(object$sma)
    if (which %in% c("routing", "both"))
        ucall["routing"] <- list(object$routing)
    eval(ucall)
}
