## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##

tf <-
    function(DATA = zoo(),
             pars = numeric(),
             delay = 0, #ihacres.getOption("delay"),
             warmup = ihacres.getOption("warmup"),
             initX = TRUE,
             na.action = na.pass,
             epsilon = TRUE)
{
    if (is.na(delay)) delay <- 0
    if (missing(DATA)) {
        ## allow to define tf without data, for later predict()ion
        DATA <- as.ts(cbind(U = rep(0, length = warmup+1),
                            Q = rep(0, length = warmup+1)))
    }
    ## get data into the right form
    DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    stopifnot(c("U","Q") %in% colnames(DATA))
    ## the model object
    obj <- list(call = match.call())
    obj$data <- DATA
    class(obj) <- c("tf")
    obj <- update(obj, pars = pars, delay = delay, warmup = warmup,
                  initX = initX, epsilon = epsilon)
    if (!inherits(obj, "tf"))
        return(obj)
    obj$call <- match.call()
    return(obj)
}

update.tf <-
    function(object, newdata = NULL, pars, ...)
{
    upcall <- match.call()
    nm <- names(upcall)
    if (!is.null(nm)) {
        nm <- nm[!(nm %in% c("", "object", "newdata"))]
        if (length(nm) > 0)
            object$call[nm] <- upcall[nm]
    }
    ## update arguments (delay, warmup, initX, epsilon)
    dots <- list(...)
    stopifnot(all(names(dots) %in%
              c("delay", "warmup", "initX", "epsilon")))
    object <- modifyList(object, dots)
    ## update parameters
    if (!missing(pars)) {
        if (length(pars) == 0) {
            object$order <- c(n=0, m=0)
            object$coefficients <- numeric()
        } else {
            pars <- tfParsConvert(pars, "a,b")
            ## check parameters
            if (isTRUE(ihacres.getOption("catch.errors"))) {
                pcheck <- try(tfParsCheck(pars),
                              silent = !ihacres.getOption("trace"))
                if (!isTRUE(pcheck)) {
                    return(pcheck)
                }
            } else {
                tfParsCheck(pars)
            }
            ## model order
            parSymbols <- gsub("_.*", "", names(pars))
            a <- pars[parSymbols == "a"]
            b <- pars[parSymbols == "b"]
            a <- stripzeros(a)
            b <- stripzeros(b, up.to = 1)
            n <- length(a)
            m <- length(b) - 1
            object$order <- c(n=n, m=m)
            extraPars <- pars[(parSymbols %in% c("a", "b")) == FALSE]
            object$coefficients <- c(a, b, extraPars)
        }
    }
    ## update DATA
    if (!is.null(newdata)) {
        newdata <- as.ts(newdata)
        stopifnot("Q" %in% colnames(newdata))
        object$data <- newdata
    }
    if (length(coef(object)) > 0) {
        ## take starting value of filter from data
        initX <- object$initX
        Xs_0 <- if (is.logical(initX)) 0 else initX
        if (isTRUE(initX)) {
            Q <- object$data[,"Q"]
            if (any(is.finite(Q[1:10]))) {
                Xs_0 <- min(Q[1:10], na.rm=TRUE)
            }
        }
        ## determine data resolution from data
        epsilon <- object$epsilon
        if (isTRUE(epsilon)) {
            epsilon <- ihacres.getOption("sim.epsilon")
            Q <- object$data[,"Q"]
            if (any(Q < epsilon, na.rm=TRUE)) {
                isNonZero <- is.finite(Q) & (Q > 0)
                if (any(isNonZero)) {
                    eps.data <- min(Q[isNonZero])
                    if (eps.data < 1)
                        epsilon <- eps.data / 10
                }
            }
            object$epsilon <- epsilon
        }
    }
    ## run model
    X <- predict(object, Xs_0 = Xs_0, epsilon = epsilon)
    object$fitted.values <- X
    object$residuals <- stripWarmup(object$data[,"Q"] - X, object$warmup)
    return(object)
}

predict.tf <-
    function(object,
             newdata = NULL,
             return_components = FALSE,
             ...)
{
    if (is.null(newdata))
        newdata <- object$data
    ## get data into the right form
    newdata <- as.ts(newdata)
    if (NCOL(newdata) > 1)
        newdata <- newdata[,"U"]
    ## simulate
    tf.sim(newdata, pars = coef(object),
                delay = object$delay,
                return_components = return_components, ...)
}

print.tf <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nUnit Hydrograph / Linear transfer function model\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    n <- x$order["n"]
    m <- x$order["m"]
    cat("Order: (n=", n, ", m=", m, ")  Delay: ", x$delay, "\n", sep="")
    cat("ARMAX Parameters:\n")
    print.default(coef(x), digits = digits, print.gap = 2)
    cat("Exponential component parameters:\n")
    pars.tv <- try(coef(x, "tau,v"), silent = TRUE)
    print.default(pars.tv, digits = digits, print.gap = 2)
    cat("TF Structure: ")
    series <- -1
    if (is.numeric(pars.tv) && ("series" %in% names(pars.tv)))
        series <- pars.tv[["series"]]
    if (n == 0) {
        cat("instantaneous")
    } else if (n == 1) {
        if (m == 0) cat("single store")
        if (m == 1) cat("single store + instantaneous in parallel")
        if (m >= 2) cat("single store + complex MA component")
    } else if (n == 2) {
        if (m == 0) cat("S * Q (two stores in series)")
        if (m == 1) cat("S + Q (two stores in parallel)")
        if (m == 2) cat("S + Q + inst. (three in parallel)")
        if (m >= 3) cat("S + Q + complex MA component")
    } else if (n == 3) {
        if (series == 0) cat("S + Q + 3 (three in parallel)")
        if (series == 1) cat("(S * 3) + Q (two in series, one in parallel)")
        if (series == 2) cat("(S + Q) * 3 (two in parallel, one in series)")
        if (series == 3) cat("S * Q * 3 in series")
        stopifnot((m == 0) == (series == 3))
    } else {
        cat("complex (high order)")
    }
    cat("\n\n")
    invisible(x)
}

coef.tf <- function(object, form = c("a,b", "tau,v", "alpha,beta"), ...)
{
    coef <- object$coefficients
    if (!missing(form)) {
        form <- match.arg(form)
        return(tfParsConvert(coef, form))
    }
    coef
}

summary.tf <-
    function(object,
             which = c("yic", "arpe", "rel.bias", "bias", "r.squared",
             "r.sq.sqrt", "r.sq.log", "r.sq.monthly",
             "ssg", "residuals"),
             ...)
{
    which <- match.arg(which, several.ok=TRUE)
    ans <- list()
    copyVars <- c("call", "order", "delay", "coefficients")
    ans[copyVars] <- object[copyVars]
    ## always want this to print with the parameter values
    ans$coef.var <- diag(vcov(object))
    ## align time series and exclude warmup period
    ## we know that object$data and fitted(object) have the same tsp
    stopifnot(tsp(object$data) == tsp(object$fitted.values))
                                        #data <- ts.intersect(obs=object$data[,"Q"], mod=fitted(object))
                                        #	warmup <- object$warmup
                                        # t_warm <- tsp(data)[1] + (warmup+1) / frequency(data)
                                        # data <- window(data, start=t_warm) # SAFE BUT SLOW
                                        #data <- data[-(1:warmup),]
                                        #data <- ts(data, t_warm, freq=frequency(data))
                                        # use vectors rather than time series, so arithmetic is faster
                                        #mod <- unclass(data[,"mod"])
                                        #	obs <- object$data[,"Q"]
                                        #	mod <- fitted(object)
                                        #	obs <- unclass(obs)[-(1:warmup)]
                                        #	mod <- unclass(mod)[-(1:warmup)]
    ## Residuals
    if ("residuals" %in% which)
        ans$residuals <- object$residuals
    ## ARPE
    if ("arpe" %in% which) {
        if (is.null(vcov(object))) ans$arpe <- NA else
        ans$arpe <- mean(ans$coef.var / (coef(object)^2))
    }
    ## YIC
    if ("yic" %in% which) {
        if (is.null(vcov(object))) ans$arpe <- NA else
        ans$arpe <- mean(ans$coef.var / (coef(object)^2))
        var.ratio <- (var(residuals(object), na.rm = TRUE) /
                      var(observed(object), na.rm = TRUE))
        ans$yic <- log(var.ratio) + log(ans$arpe)
    }
    ## Steady State Gain
    if ("ssg" %in% which) {
        ans$ssg <- ssg.tf.coef(coef(object))
    }
    ## call perfStats for the rest
    DATA <- ts.intersect(Q=object$data[,"Q"], U=object$data[,"U"], X=fitted(object))
    ans <- c(ans, as.list(perfStats(DATA, warmup=object$warmup, which=which, ...)))
    class(ans) <- "summary.tf"
    ans
}

print.summary.tf <-
    function(x, digits = max(3, getOption("digits") - 3), ...,
             print.call = TRUE)
{
    with(x, {
        if (print.call) {
            cat("\nCall:\n")
            print(call)
        }
        cat("\n")
        n <- order["n"]
        m <- order["m"]
        cat("Order: (n=", n, ", m=", m, ")  Delay: ", delay, "\n", sep="")
        cat("ARMAX Parameters:\n")
        theta.etc <- coefficients
        if (length(coef.var) > 0) {
            theta.etc <- rbind(theta.etc, std.dev.=sqrt(coef.var),
                               deparse.level=0)
        }
        print.default(theta.etc, digits=digits, print.gap=2)
        cat("\n")
        if (!is.null(x$residuals)) {
            cat("Residuals (obs - mod): \n")
            print(summary(residuals))
        }
        if (!is.null(x$yic))
            cat("YIC:", format(yic, digits=digits), "\n")
        if (!is.null(x$arpe))
            cat("ARPE (%):", format(arpe * 100, digits=digits), "\n")
        if (!is.null(x$r.squared))
            cat("R Squared:", format(r.squared, digits=digits), "\n")
        if (!is.null(x$r.sq.sqrt))
            cat("R Squared sqrt:", format(r.sq.sqrt, digits=digits), "\n")
        if (!is.null(x$r.sq.log))
            cat("R Squared log:", format(r.sq.log, digits=digits), "\n")
        if (!is.null(x$r.sq.monthly))
            cat("R Squared monthly:", format(r.sq.monthly, digits=digits), "\n")
        if (!is.null(x$bias))
            cat("Bias (units, total):", format(bias, digits=digits), "\n")
        if (!is.null(x$rel.bias))
            cat("Rel. Bias (%):", format(rel.bias * 100, digits=digits), "\n")
        if (!is.null(x$U1))
            cat("Error corr. with lag input (U1):", format(U1, digits=digits), "\n")
        if (!is.null(x$X1))
            cat("Error corr. with lag output (X1):", format(X1, digits=digits), "\n")
        if (!is.null(x$ssg))
            cat("Steady State Gain:", format(ssg, digits=digits), "\n")
        if (!is.null(x$obsgain))
            cat("Observed Gain (Q/U):", format(obsgain, digits=digits), "\n")
        cat("\n")
    })
    invisible(x)
}

ssg.tf <-
    function(object, ...)
{
    theta <- coef(object)
#    if (length(theta) == 0) {
#        ## estimate scale by mass balance
#        Q <- object$data[,"Q"]
#        U <- object$data[,"U"]
#        ok <- complete.cases(Q, U)
#        return(sum(U[ok]) / sum(Q[ok]))
#    } else {
        ssg.tf.coef(theta)
#    }
}


normalise.tf <-
    function(object, ...)
{
    ## TODO: do not run model again, just scale output?
    update(object, pars = normalise.tf.coef(coef(object)))
}

fitted.tf <- function(object, ..., all = FALSE)
{
    if (all)
        return(object$fitted.values)
    stripWarmup(object$fitted.values, object$warmup)
}

residuals.tf <- function(object, ..., all = FALSE)
{
    tmp <- object$data[,"Q"] - object$fitted.values
    if (all) return(tmp)
    stripWarmup(tmp, object$warmup)
}

observed.tf <- function(object, ..., all = FALSE)
{
    ## observed.default will work, but this may be slightly faster
    if (all)
        return(object$data[,"Q"])
    stripWarmup(object$data[,"Q"], object$warmup)
}

getU.tf <- function(object, ..., all = FALSE)
{
    if (all)
        return(object$data[,"U"])
    stripWarmup(object$data[,"U"], object$warmup)
}

vcov.tf <- function(object, ...)
    object$cov.mat ## may be NULL

deviance.tf <- stats:::deviance.lm

isValidModel <- function(object, ...)
    UseMethod("isValidModel")

isValidModel.default <- function(object, ...)
    return(FALSE)

isValidModel.tf <- function(object, ...)
    (inherits(object, "tf") &&
     (isTRUE(try(tfParsCheck(coef(object)), silent = TRUE))))
