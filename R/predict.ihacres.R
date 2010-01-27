## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


predict.ihacres <-
    function(object, newdata = NULL,
             with_routing = TRUE,
             return_state = FALSE,
             return_components = FALSE, ...)
{
    if (is.null(newdata)) {
        newdata <- object$data
    } else {
        newdata <- as.ts(newdata)
    }
    ## check that parameters are fully specified
    which <- if (with_routing) "both" else "sma"
    if (!isFullySpecified(object, which = which))
        stop("model parameters are not fully specified")
    ## construct call to simulation function
    simpars <- coef(object, which = which)
    ucall <- quote(ihacres.sim(newdata))
    ucall <- as.call(c(as.list(ucall),
                       as.list(simpars)))
    ucall["sma"] <- list(object$sma)
    if (with_routing) {
        ucall["routing"] <- list(object$routing)
    } else {
        ucall["routing"] <- list(NULL)
    }
    if (return_state)
        ucall$return_state <- TRUE
    eval(ucall)
}
