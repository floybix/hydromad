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
    dots <- list(...)
    if (!is.null(dots$warmup)) stop("'warmup' can not be given here")
    if (!is.null(dots$normalise)) stop("'normalise' can not be given here")
    model <- tf.ls.fit(DATA, order = order, delay = delay,
                       warmup = 0, normalise = FALSE, ...)
    if (!inherits(model, "tf"))
        return(model)
    n <- order[1]
    poles <- arToPoles(coef(model)[seq_len(n)])
    eps <- sqrt(.Machine$double.eps)
    badpoles <- (Re(poles) < -eps) | (abs(Im(poles)) > eps)
    if (any(badpoles)) {
        model <- tfFitWithPoleConstraints(DATA, tf.fit = tf.ls.fit, poles = poles,
                               order = order, delay = delay,
                               warmup = 0, normalise = FALSE, ...,
                               method = method, control = control, hessian = hessian)
    }
    if (!isValidModel(model))
        return(model)

    model$coefficients <- c(coef(model, "tau,v"),
                            delay = model$delay)
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
    dots <- list(...)
    if (!is.null(dots$warmup)) stop("'warmup' can not be given here")
    if (!is.null(dots$normalise)) stop("'normalise' can not be given here")
    model <- tf.sriv.fit(DATA, order = order, delay = delay,
                         warmup = 0, normalise = FALSE, ...)
    if (!inherits(model, "tf"))
        return(model)
    n <- order[1]
    poles <- arToPoles(coef(model)[seq_len(n)])
    eps <- sqrt(.Machine$double.eps)
    badpoles <- (Re(poles) < -eps) | (abs(Im(poles)) > eps)
    if (any(badpoles)) {
        model <- tfFitWithPoleConstraints(DATA, tf.fit = tf.sriv.fit, poles = poles,
                               order = order, delay = delay,
                               warmup = 0, normalise = FALSE, ...,
                               method = method, control = control, hessian = hessian)
    }
    if (!isValidModel(model))
        return(model)
    #model$vcov <- vcov(model)
    model$coefficients <- c(coef(model, "tau,v"),
                            delay = model$delay)
    model
}


expuh.inverse.fit <-
    function(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             ...)
{
    model <- tf.inverse.fit(DATA, order = order, delay = delay,
                            ...)
    if (!inherits(model, "tf"))
        return(model)
    #model$vcov <- vcov(model)
    model$coefficients <- c(coef(model, "tau,v"),
                            delay = model$delay)
    model$fitted.values <- NULL
    model$residuals <- NULL
    model
}

tfFitWithPoleConstraints <- function(DATA, tf.fit, poles, ...,
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
            thisMod <- try(tf.fit(DATA, ..., fixed.ar = polesToAr(pol)))
        } else {
            thisMod <- tf.fit(DATA, ..., fixed.ar = polesToAr(pol))
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
