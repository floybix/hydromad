## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


predict.hydromad <-
    function(object, newdata = NULL,
             which = c("both", "sma", "routing"),
             ..., all = TRUE,
             return_state = FALSE,
             return_components = FALSE)
{
    which <- match.arg(which)
    if (is.null(newdata)) {
        if (which == "routing") {
            newdata <- object$U
        } else {
            newdata <- object$data
        }
    } else {
        newdata <- as.zooreg(newdata)
    }
    DATA <- newdata
    ## check that parameters are fully specified
    if (!isFullySpecified(object, which = which))
        stop("model parameters are not fully specified")
    
    ## construct calls to simulation functions
    sma <- object$sma
    routing <- object$routing
    if (which == "routing") sma <- NULL
    if (which == "sma") routing <- NULL
    
    if (is.character(sma)) {
        ## construct call to SMA simulation function
        sma.fun <- paste(sma, ".sim", sep = "")
        sma.args <- coef(object, which = "sma", etc = TRUE)
        ucall <- as.call(c(list(as.symbol(sma.fun),
                                quote(DATA)),
                           sma.args))
        if (return_state)
            ucall$return_state <- TRUE
        ## calculate U
        U <- eval(ucall)
        if (return_state) {
            S <- U
            if (NCOL(S) > 1) {
                stopifnot("U" %in% colnames(S))
                U <- S[,"U"]
            }
        }
    } else {
        ## take observed rainfall P as effective rainfall U
        if (NCOL(DATA) > 1) {
            stopifnot("P" %in% colnames(DATA))
            U <- DATA[,"P"]
        } else {
            U <- DATA
        }
        if (return_state)
            S <- U
    }
    ## handle routing
    if (is.character(routing)) {
        ## construct call to routing simulation function
        r.fun <- paste(routing, ".sim", sep = "")
        r.args <- coef(object, which = "routing", etc = TRUE)
        ## TODO: take starting value of filter from data?
        rcall <- as.call(c(list(as.symbol(r.fun),
                            quote(U)),
                           r.args))
        if (return_components)
            rcall$return_components <- TRUE
        Q <- eval(rcall)
    } else {
        ## no routing
        Q <- U
    }
    if (return_state) {
        if (is.null(routing)) {
            ans <- S
        } else {
            ans <- cbind(S, Q)
            if (length(colnames(S)) > 0)
                colnames(ans)[1:NCOL(S)] <- colnames(S)
        }
    } else {
        ans <- Q
    }
    if (all) ans else stripWarmup(ans, object$warmup)
}
