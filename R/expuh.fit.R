## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

expuh.ls.fit <-
    function(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             ...,
             method = hydromad.getOption("optim.method"),
             control = hydromad.getOption("optim.control"),
             hessian = TRUE)
{
    model <- armax.ls.fit(DATA, order = order, delay = delay, ...)
    if (!inherits(model, "tf"))
        return(model)
    n <- order[1]
    poles <- arToPoles(coef(model)[seq_len(n)])
    eps <- sqrt(.Machine$double.eps)
    badpoles <- (Re(poles) < -eps) | (abs(Im(poles)) > eps)
    if (any(badpoles)) {
        model <-
            fitWithPoleConstraints(DATA, fitfun = armax.ls.fit, poles = poles,
                                   order = order, delay = delay, ...,
                                   method = method, control = control, hessian = hessian)
    }
    if (!inherits(model, "tf"))
        return(model)
    model$coefficients <- coef(model, "tau,v")
    model
}


expuh.sriv.fit <-
    function(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             ...,
             method = hydromad.getOption("optim.method"),
             control = hydromad.getOption("optim.control"),
             hessian = TRUE)
{
    model <- armax.sriv.fit(DATA, order = order, delay = delay, ...)
    if (!inherits(model, "tf"))
        return(model)
    n <- order[1]
    poles <- arToPoles(coef(model)[seq_len(n)])
    eps <- sqrt(.Machine$double.eps)
    badpoles <- (Re(poles) < -eps) | (abs(Im(poles)) > eps)
    if (any(badpoles)) {
        model <-
            fitWithPoleConstraints(DATA, fitfun = armax.sriv.fit, poles = poles,
                                   order = order, delay = delay, ...,
                                   method = method, control = control, hessian = hessian)
    }
    if (!inherits(model, "tf"))
        return(model)
    model$coefficients <- coef(model, "tau,v")
    model
}


expuh.inverse.fit <-
    function(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             ...)
{
    ## TODO: can do this directly?
    model <- armax.inverse.fit(DATA, order = order, delay = delay, ...)
    if (!inherits(model, "tf"))
        return(model)
    #model$vcov <- vcov(model)
    model$coefficients <- coef(model, "tau,v")
    model
}

fitWithPoleConstraints <- function(DATA, fitfun, poles, ...,
                                   method, control, hessian)
{
    ## non-physical roots; constrain them and fit again
    poles <- abs(Re(poles))
    poles <- pmax(poles, 1e-5)
    logpol0 <- log(poles)
    model <- structure("failed to fit expuh routing with constrained roots",
                       class = "try-error")
    bestVal <- Inf
    optFun <- function(logpol) {
        pol <- exp(logpol)
        if (isTRUE(hydromad.getOption("catch.errors"))) {
            thisMod <- try(fitfun(DATA, ..., fixed.ar = polesToAr(pol)))
        } else {
            thisMod <- fitfun(DATA, ..., fixed.ar = polesToAr(pol))
        }
        if (!isValidModel(thisMod))
            return(NA)
        val <- objFunVal(thisMod)
        if (val < bestVal)
            model <<- thisMod
        val
    }
    if (!isTRUE(hydromad.getOption("catch.errors.optim")))
        try <- force ## i.e. skip the try()
    ans <- try(optim(logpol0, optFun, method = method,
                     control = control, hessian = hessian))
    if (inherits(ans, "try-error")) {
        return(ans)
    }
    if (hessian) {} ## TODO
    if (ans$convergence != 0) {
        msg <- if (ans$convergence == 1) {
            "optim() reached maximum iterations"
        } else {
            paste("optim() returned convergence code",
                  ans$convergence)
        }
        if (!isTRUE(hydromad.getOption("quiet"))) {
            warning(msg)
        }
    }
    model
}
