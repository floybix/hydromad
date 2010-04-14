## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

coef.hydromad <-
    function(object, which = c("both", "sma", "routing"), ...,
             warn = TRUE, etc = FALSE)
{
    which <- match.arg(which)
    parlist <- object$parlist
    if (etc == FALSE) {
        ## only return (numeric) parameters, ignore other objects
        ok <- unlist(lapply(parlist, function(x)
                 {
                     (!inherits(x, "AsIs")) &&
                     (is.numeric(x) ||
                      (is.atomic(x) && all(is.na(x))))
                 }))
        parlist <- parlist[ok]
    }
    if (which == "both") {
        if (etc) {
            return(parlist)
        }
        if (any(sapply(parlist, length) > 1)) {
            if (warn)
                warning("parameters not fully specified, returning list")
            return(parlist)
        } else {
            return(unlist(parlist))
        }
    }
    ## work out which arguments go to SMA function
    sma.argnames <- names(object$sma.formals)
    forSMA <- names(parlist) %in% sma.argnames
    ## work out which arguments go to routing function
    routing <- object$routing
    r.argnames <- NULL
    if (is.character(routing)) {
        r.fun <- paste(routing, ".sim", sep = "")
        r.argnames <- names(formals(r.fun))
    }
    forRouting <- names(parlist) %in% r.argnames
    ## resolve ambiguities / arguments to be passed through
    unmatched <- (!forSMA) & (!forRouting)
    unmatchedOK <- FALSE
    if ("..." %in% sma.argnames) {
        forSMA <- forSMA | unmatched
        unmatchedOK <- TRUE
    }
    if ("..." %in% r.argnames) {
        forRouting <- forRouting | unmatched
        unmatchedOK <- TRUE
    }
    if (any(unmatched) && !unmatchedOK)
        warning("unrecognised parameters: ",
                toString(names(parlist)[unmatched]))
    if (which == "sma")
        result <- parlist[forSMA]
    if (which == "routing")
        result <- parlist[forRouting]
    if (etc) {
        return(result)
    }
    if (any(sapply(result, length) > 1)) {
        if (warn)
            warning("parameters not fully specified, returning list")
        return(result)
    } else {
        return(unlist(result))
    }
}

"coef<-" <- function(object, ..., value)
    UseMethod("coef<-")

"coef<-.hydromad" <-
    function(object, ..., value)
{
    pars <- as.list(value)
    ## all elements must have names
    if (length(names(pars)) < length(pars) ||
        any(names(pars) == ""))
    {
        stop("elements to assign to coef() must be named")
    }
    ## find existing parameters not named in 'value', to remove:
    curNames <- names(coef(object, ..., warn = FALSE))
    remNames <- curNames[curNames %in% names(pars) == FALSE]
    ## set them to NULL to delete them from the coef() list:
    remList <- rep(list(NULL), length(remNames))
    names(remList) <- remNames
    pars <- c(pars, remList)
    ## use 'update' to make the changes
    do.call("update", c(quote(object), pars))
}
