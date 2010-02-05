## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


uh.sim <-
    function(U,
             a_1 = 0, a_2 = 0, a_3 = 0, a_4 = 0, a_5 = 0, a_6 = 0, a_7 = 0, a_8 = 0,
             b_0 = 1, b_1 = 0, b_2 = 0, b_3 = 0, b_4 = 0, b_5 = 0, b_6 = 0, b_7 = 0,
             delay = 0, X_0 = 0)
{
    pars <- unlist(data.frame(a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8,
                              b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7))
    tf.sim(U, pars = pars, delay = delay, Xs_0 = X_0)
}


ssg.uh <- function(theta)
{
    ssg.tf.coef(theta)
}

normalise.uh <- function(theta)
    normalise.tf.coef(theta)


uh.ls.fit <-
    function(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             ...)
{
    dots <- list(...)
    if (!is.null(dots$warmup)) stop("'warmup' can not be given here")
    if (!is.null(dots$normalise)) stop("'normalise' can not be given here")
    model <- tf.ls.fit(DATA, order = order, delay = delay,
              warmup = 0, normalise = FALSE, ...)
    if (!inherits(model, "tf"))
        return(model)
    model$coefficients <- c(coef(model),
                            delay = model$delay)
    model
}

uh.sriv.fit <-
    function(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             ...)
{
    dots <- list(...)
    if (!is.null(dots$warmup)) stop("'warmup' can not be given here")
    if (!is.null(dots$normalise)) stop("'normalise' can not be given here")
    model <- tf.sriv.fit(DATA, order = order, delay = delay,
                warmup = 0, normalise = FALSE, ...)
    if (!inherits(model, "tf"))
        return(model)
    model$coefficients <- c(coef(model),
                            delay = model$delay)
    model
}

uh.inverse.fit <-
    function(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             ...)
{
    model <- tf.inverse.fit(DATA, order = order, delay = delay,
                            ...)
    if (!inherits(model, "tf"))
        return(model)
    model$fitted.values <- NULL
    model$residuals <- NULL
    model
}
