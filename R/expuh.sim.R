## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


expuh.sim <-
    function(U,
             tau_s = 0, tau_q = 0, tau_3 = 0,
             v_s = 1, v_3 = 0, v_q = NA,
             series = 0, lambda = 0, loss = 0,
             delay = 0, Xs_0 = 0, Xq_0 = 0,
             return_components = FALSE)
{
    delay <- round(delay)
    if (is.na(v_q))
        v_q <- max(0, min(1, 1 - v_s - v_3))
    pars <- c(tau_s = tau_s,
              tau_q = tau_q,
              tau_3 = tau_3,
              v_s = v_s,
              v_q = v_q,
              v_3 = v_3,
              if (series != 0) c(series = series),
              if (lambda != 0) c(lambda = lambda),
              if (loss != 0) c(loss = loss))
    tf.sim(U, pars = pars, delay = delay,
           Xs_0 = Xs_0, Xq_0 = Xq_0,
           return_components = return_components)
}


ssg.expuh <- function(theta)
{
    if (length(theta) == 0)
        return(1)
    theta <- tfParsConvert(theta, "a,b")
    ssg.tf.coef(theta)
}

normalise.expuh <- function(theta)
{
#    v_s <- 1
#    v_3 <- 0
#    v_q <- NA
#    if ("v_s" %in% names(theta))
#        v_s <- theta[["v_s"]]
#    if ("v_3" %in% names(theta))
#        v_3 <- theta[["v_3"]]
#    if ("v_q" %in% names(theta))
#        v_q <- theta[["v_q"]]
#    if (is.na(v_q)) {
#        ## unit volume is enforced anyway
#        return(theta)
#    }
#    vv <- v_s + v_q + v_3

    theta <- tfParsConvert(theta, "a,b")
    tmp <- normalise.tf.coef(theta)
    tfParsConvert(tmp, "tau,v")
}

expuh.ls.fit <-
    function(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             warmup = stop("ignored"),
             normalise = stop("ignored"),
             ...,
             method = hydromad.getOption("optim.method"),
             control = hydromad.getOption("optim.control"),
             hessian = TRUE)
{
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
             warmup = stop("ignored"),
             normalise = stop("ignored"),
             ...,
             method = hydromad.getOption("optim.method"),
             control = hydromad.getOption("optim.control"),
             hessian = TRUE)
{
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
