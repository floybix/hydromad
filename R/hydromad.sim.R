## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


hydromad.sim <-
    function(DATA,
             ...,
             sma = hydromad.getOption("sma"),
             routing = hydromad.getOption("routing"),
             return_state = FALSE,
             return_components = FALSE)
{
    DATA <- as.ts(DATA)
    ## dots `...` may contain arguments for sma and/or routing
    dots <- list(...)
    ## work out which arguments go to SMA function
    sma.argnames <- NULL
    if (is.character(sma)) {
        sma.fun <- paste(sma, ".sim", sep = "")
        sma.argnames <- names(formals(sma.fun))
    }
    forSMA <- names(dots) %in% sma.argnames
    ## work out which arguments go to routing function
    r.argnames <- NULL
    if (is.character(routing)) {
        r.fun <- paste(routing, ".sim", sep = "")
        r.argnames <- names(formals(r.fun))
    }
    forRouting <- names(dots) %in% r.argnames
    ## resolve ambiguities / arguments to be passed through
    unmatched <- (!forSMA) & (!forRouting)
    if ("..." %in% sma.argnames) {
        forSMA <- forSMA | unmatched
    }
    if ("..." %in% r.argnames) {
        forRouting <- forRouting | unmatched
    }
    unmatched <- (!forSMA) & (!forRouting)
    if (any(unmatched)) {
        stop("unused arguments: ", toString(names(dots)[unmatched]))
    }
    if (is.character(sma)) {
        ## construct call to SMA simulation function
        ucall <- as.call(c(list(as.symbol(sma.fun),
                                quote(DATA)),
                           dots[forSMA]))
        if (return_state)
            ucall$return_state <- TRUE
        ## calculate U
        U <- eval(ucall)
        if (return_state) {
            if (is.null(routing))
                return(U)
            S <- U
            U <- S[,"U"]
        }
    } else if (is.null(sma)) {
        ## take observed rainfall P as effective rainfall U
        if (NCOL(DATA) > 1) {
            U <- DATA[,"P"]
        } else {
            U <- DATA
        }
        if (return_state)
            S <- U
    } else {
        stop("unrecognised value of 'sma'")
    }
    ## handle routing
    if (is.character(routing)) {
        ## construct call to routing simulation function

        ## TODO: take starting value of filter from data

        rcall <- as.call(c(list(as.symbol(r.fun),
                            quote(U)),
                       dots[forRouting]))
        if (return_components)
            rcall$return_components <- TRUE
        Q <- eval(rcall)
    } else if (is.null(routing)) {
        ## no routing
        Q <- U
    } else {
        stop("unrecognised value of 'routing'")
    }
    if (return_state) {
        ans <- cbind(S, Q)
        if (length(colnames(S)) > 0)
            colnames(ans)[1:NCOL(S)] <- colnames(S)
        return(ans)
    } else {
        return(Q)
    }
}

