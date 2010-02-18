## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

hydromad <-
    function(DATA = zoo(),
             ...,
             sma = hydromad.getOption("sma"),
             routing = hydromad.getOption("routing"),
             rfit = NULL,
             weights = NULL,
             warmup = hydromad.getOption("warmup"))
{
    ## TODO: keep zoo attributes?
    DATA <- as.ts(DATA)
    ## dots `...` may contain arguments for sma and/or routing.
    ## parlist stores these -- and they may be ranges.
    ## take defaults from hydromad.options()
    parlist <- list()
    if (is.character(routing)) {
        parlist <-
            modifyList(parlist, as.list(hydromad.getOption(routing)))
    }
    if (is.character(sma)) {
        parlist <-
            modifyList(parlist, as.list(hydromad.getOption(sma)))
    }
    parlist <- modifyList(parlist, list(...))
    ## create the model object
    obj <- list(call = match.call())
    class(obj) <- unique(c(sma, "hydromad"))
    obj$sma <- sma
    if (is.character(sma)) {
        obj$sma.fun <- paste(sma, ".sim", sep = "")
        force(get(obj$sma.fun, mode = "function"))
        obj$sma.args <- formals(obj$sma.fun)
    }
    obj$parlist <- parlist
    obj <- update(obj, newdata = DATA,
                  routing = routing, rfit = rfit,
                  warmup = warmup, weights = weights)
    obj$call <- match.call() ## reset call after update()
    return(obj)
}

isFullySpecified <- function(object, ...)
    !is.list(coef(object, ..., warn = FALSE))

coef.hydromad <-
    function(object, which = c("both", "sma", "routing"), ..., warn = TRUE)
{
    which <- match.arg(which)
    parlist <- object$parlist
    if (which == "both") {
        if (any(sapply(parlist, length) > 1)) {
            if (warn)
                warning("parameters not fully specified, returning list")
            return(parlist)
        } else {
            return(unlist(parlist))
        }
    }
    ## work out which arguments go to SMA function
    sma.argnames <- names(object$sma.args)
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
    if (any(sapply(result, length) > 1)) {
        if (warn)
            warning("parameters not fully specified, returning list")
        return(result)
    } else {
        return(unlist(result))
    }
}

fitted.hydromad <- function(object, ..., U = FALSE, all = FALSE)
{
    if (is.null(object$routing))
        U <- TRUE
    tmp <- if (U) object$U else object$fitted.values
    if (length(tmp) == 0)
        return(tmp)
    if (all) return(tmp)
    stripWarmup(tmp, object$warmup)
}

residuals.hydromad <- function(object, ..., all = FALSE)
{
    f <- fitted(object, all = TRUE)
    if (length(f) == 0) return(f)
    tmp <- (object$data[,"Q"] - f)
    if (all) return(tmp)
    stripWarmup(tmp, object$warmup)
}

observed.hydromad <- function(object, ..., item = "Q", all = FALSE)
{
    ## observed.default will work (for Q), but this may be slightly faster
    tmp <- object$data[,item]
    if (all) return(tmp)
    stripWarmup(tmp, object$warmup)
}

vcov.hydromad <- function(object, ...)
{
    cov.mat <- object$cov.mat ## may be NULL
    rcov <- object$vcov.rfit
    ans <- cov.mat
    if (!is.null(rcov)) {
        if (is.null(ans)) {
            ans <- rcov
        } else {
            ## merge the two matrices
            tmp.right <- rbind(matrix(ncol = ncol(rcov), nrow = nrow(cov.mat)),
                               rcov)
            ans <- rbind(cov.mat, matrix(ncol = ncol(cov.mat), nrow = nrow(rcov)))
            ans <- cbind(ans, tmp.right)
            rownames(ans) <- colnames(ans) <-
                c(rownames(cov.mat), rownames(rcov))
        }
    }
    ans
}

deviance.hydromad <- stats:::deviance.lm

print.hydromad <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    rx <- x$data
    routingStr <- if (is.null(x$routing)) "NULL" else toString(x$routing)
    cat(paste("\nHydromad model with SMA class \"", class(x)[1],
              "\" and routing \"", routingStr, "\":\n", sep = ""))
    cat(paste("Start = ", index2char(index(rx)[1], frequency(rx)),
              ", End = ", index2char(index(rx)[NROW(rx)], frequency(rx)),
              "\n", sep = ""))
    cat("Last updated: ", format(x$last.updated), "\n")
    cat("\nCall:\n")
    print(x$call)
    if (!is.null(x$objective)) {
        cat("Objective function:\n")
        print(x$objective)
    }
    cat("\n")
    for (which in c("sma", "routing")) {
        if (!is.null(x[[which]])) {
            if (which == "sma") {
                cat("SMA Parameters:\n")
            } else {
                cat("Routing Parameters:\n")
            }
            if (isFullySpecified(x, which = which)) {
                ## all unique parameter values
                coefx <- coef(x, which = which)
                print(coefx, digits = digits, quote = FALSE, print.gap = 2)
                if (which == "routing") {
                    tmp <- describeTF(coefx)
                    if (!is.null(tmp) && !is.na(tmp))
                        cat("TF Structure: ", tmp, "\n")
                }
            } else {
                ## one or more parameters specified as ranges only
                parlist <- coef(x, which = which, warn = FALSE)
                print(data.frame(lower = sapply(parlist, min),
                                 upper = sapply(parlist, max),
                                 ` ` = sapply(parlist, function(p)
                             {
                                 ## mark fixed parameters
                                 if (isTRUE(diff(range(p)) == 0)) '(==)' else ''
                             }), check.names = FALSE),
                      digits = digits)
            }
            cat("\n")
        }
    }
    if (!is.null(x$rfit)) {
        cat("Routing fit spec.: ",
            toString(deparse(x$rfit, control = c(), width = 500),
                     width = 60), "\n")
    }
    if (!is.null(x$msg)) {
        cat("\nMessage: ", toString(x$msg), "\n")
    }
    invisible(x)
}

isValidModel <- function(object, ...)
    UseMethod("isValidModel")

isValidModel.default <- function(object, ...)
    return(FALSE)

isValidModel.hydromad <- function(object, ...)
    is.numeric(fitted(object, all = TRUE))
