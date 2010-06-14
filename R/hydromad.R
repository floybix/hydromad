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
    ## create the model object
    obj <- list(call = match.call())
    class(obj) <- "hydromad"
    ## dots `...` may contain arguments for sma and/or routing.
    ## parlist stores these -- and they may be ranges.
    ## update() takes default parameter ranges/values from hydromad.options().
    parlist <- list(...)
    obj$parlist <- parlist
    obj <- update(obj, newdata = DATA, sma = sma,
                  routing = routing, rfit = rfit,
                  warmup = warmup, weights = weights)
    obj$call <- match.call() ## reset call after update()
    return(obj)
}

isFullySpecified <- function(object, ...)
    !is.list(coef(object, ..., warn = FALSE))

fitted.hydromad <- function(object, ..., U = FALSE, all = FALSE)
{
    if (is.null(object$routing))
        U <- TRUE
    tmp <- if (U) object$U else object$fitted.values
    if (length(tmp) == 0)
        return(tmp)
    if (all) tmp else stripWarmup(tmp, object$warmup)
}

residuals.hydromad <- function(object, ..., all = FALSE)
{
    f <- fitted(object, all = TRUE)
    if (length(f) == 0) return(f)
    tmp <- (object$data[,"Q"] - f)
    if (all) tmp else stripWarmup(tmp, object$warmup)
}

observed.hydromad <- function(object, ..., item = "Q", all = FALSE)
{
    ## observed.default will work (for Q), but this may be slightly faster
    if (is.character(item))
        if (!all(item %in% colnames(object$data)))
            return(NULL)
    tmp <- object$data[,item]
    if (all) tmp else stripWarmup(tmp, object$warmup)
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

logLik.hydromad <-
    function(object, loglik = hydromad.getOption("loglik"), ...)
{
    val <- objFunVal(object, objective = loglik, ...)
    ## TODO: for a fitted model we do not know how many parameters were fitted
    ## guess:
    attr(val, "df") <- length(coef(object))
    attr(val, "nobs") <- attr(val, "nall") <- length(fitted(object, all = TRUE))
    class(val) <- "logLik"
    val
}

deviance.hydromad <- stats:::deviance.lm

print.hydromad <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat(paste("\n",
              'Hydromad model with "', format(x$sma), '" SMA',
              ' and "', format(x$routing), '" routing:', "\n", sep = ''))
    rx <- x$data
    cat(paste("Start = ", index2char(index(rx)[1], frequency(rx)),
              ", End = ", index2char(index(rx)[NROW(rx)], frequency(rx)),
              "\n", sep = ""))
    cat("Last updated: ", format(x$last.updated), "\n")
    cat("\nCall:\n")
    print(x$call)
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
            #cat("\n")
        }
    }
    if (!is.null(x$rfit)) {
        cat("Routing fit spec.: ",
            toString(deparse(x$rfit, control = c(), width = 500),
                     width = getOption("width")), "\n")
    }
    if (!is.null(x$msg)) {
        cat("\nMessage: ", toString(x$msg), "\n")
    }
    if (!is.null(x$fit.call)) {
        cat("\nFit: (see $fit.result)\n")
        print(x$fit.call)
        cat("Function evaluations: ", x$funevals, " in ",
            x$timing[3], " seconds", "\n")
        cat("Objective: ",
            toString(paste(deparse(x$objective, control = c(),
                             width = 60), collapse = "\n"),
                     width = getOption("width")))
   }
    if (length(x$info.rfit) > 0) {
        cat("\nRouting fit info: ",
            toString(deparse(x$info.rfit, control = c(), width = 500),
                     width = getOption("width")), "\n")
    }
    invisible(x)
}

str.hydromad.runlist <- 
    function(object, ...)
{
    cat("\nList of Hydromad model runs:\n")
    str(lapply(object, function(obj) {
        if (!is.null(obj$msg)) {
            list(call = obj$call, message = obj$msg)
        } else {
            obj$call
        }
    }))
    invisible()
}

isValidModel <- function(object, ...)
    UseMethod("isValidModel")

isValidModel.default <- function(object, ...)
    return(FALSE)

isValidModel.hydromad <- function(object, ...)
    is.numeric(fitted(object, all = TRUE))
