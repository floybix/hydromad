## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

tf.lambda.fit <-
    function(DATA = list(U=, Q=),
             order = c(m = 2, n = 1), ## ignored
             delay = hydromad.getOption("delay"),
             which.pars = c("lambda", "tau_s", "tau_q"), ## tau_q or v_s?
             lambda.init = -0.25,
             fit.method = hydromad.getOption("inverse.fit.method"),
             ...,
             uhParTrans = alist(v_s = list(trans = log,
                                inverse = function(x) min(exp(x), 5))),
             objective = hydromad.getOption("objective"),
             optim.method = hydromad.getOption("optim.method"),
             optim.control = hydromad.getOption("optim.control"),
             trace = hydromad.getOption("trace"),
             na.action = na.pass,
             hessian = FALSE)
{
    if (inherits(objective, "formula"))
        objective <- objective[[2]]
    stopifnot(is.language(objective))

    ## TODO: merge this function with hydromad.fit.default?

    optim.control <- modifyList(hydromad.getOption("optim.control"),
                                optim.control)
    ## get data into the right form
    if (is.list(DATA)) DATA <- do.call(ts.intersect, lapply(DATA, as.ts))
    if (!is.ts(DATA)) DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    stopifnot(c("U","Q") %in% colnames(DATA))
    if (is.na(delay))
        delay <- estimateDelay(DATA, plot=FALSE)

    init.mod <- tf.fit(DATA,
                           order = c(n=2, m=1), delay = delay,
                           method = fit.method,
                           ...)
    ## initial parameter estimates
    v_s_0 <- 1.0
    init.pars <- coef(init.mod, "tau,v")[c("tau_s", "tau_q")]
    ## it is possible that tau_q here is NA
    names(init.pars) <- c("tau_s", "tau_q")
    init.pars[is.na(init.pars)] <- 0
    init.pars["v_s"] <- v_s_0
    init.pars["lambda"] <- lambda.init

    flowpars <- init.pars[which.pars]

    ## set up uh parameter transformations
    uhParTrans <- modifyList(lapply(hydromad.getOption("uhParTrans"), eval),
                             lapply(uhParTrans, eval))
    ## apply transformations
    for (nm in names(flowpars)) {
        ## check for exact match first, then partial match
        fn <- uhParTrans[[nm]]$trans
        if (is.null(fn))
            fn <- uhParTrans[[sub("_.*", "", nm)]]$trans
        if (is.function(fn))
            flowpars[[nm]] <- fn(flowpars[[nm]])
    }
    ## some transformations can induce Inf or -Inf, e.g. log(0)
    make.finite <- function(x)
        if (is.infinite(x)) sign(x) * 100 else x
    flowpars <- lapply(flowpars, make.finite)

    generateModel <- function(pars) {
        flowpp <- as.list(pars)
        ## apply inverse-transform function to each parameter
        for (nm in names(flowpp)) {
            ## check for exact match first, then partial match
            fn <- uhParTrans[[nm]]$inverse
            if (is.null(fn))
                fn <- uhParTrans[[sub("_.*", "", nm)]]$inverse
            if (is.function(fn))
                flowpp[[nm]] <- fn(flowpp[[nm]])
        }
        flowpp <- modifyList(as.list(init.pars), flowpp)
        if (flowpp$tau_q > flowpp$tau_s)
            flowpp[c("tau_q", "tau_s")] <-
                rev(flowpp[c("tau_q", "tau_s")])
        update(init.mod, pars = flowpp)
    }
    xf_optim <- function(pars) {
        #mod <- try(generateModel(pars), silent=TRUE)
        mod <- generateModel(pars)
        if (!isValidModel(mod)) {
            #warning(mod)
            return(NA)
        }
        ## these can be referred to in `objective`
        Q <- observed(mod)
        X <- fitted(mod)
        U <- getU(mod)
        eval(objective)
    }
    ans <- optim(unlist(flowpars), xf_optim, method = optim.method,
                 control = optim.control, hessian = hessian)
    if (hessian) return(ans)
    ## return tf object
    mod <- generateModel(ans$par)
    mod$call <- match.call()
    mod
}


tf.lambda.part.fit <-
    function(DATA = list(U=, Q=),
             order = c(m = 2, n = 1), ## ignored
             delay = hydromad.getOption("delay"),
             which.pars = c("lambda", "tau_s", "tau_q"), ## tau_q or v_s?
             lambda.init = -0.25,
             fit.method = hydromad.getOption("inverse.fit.method"),
             ...,
             uhParTrans = alist(v_s = list(trans = log,
                                inverse = function(x) min(exp(x), 5))),
             objective = hydromad.getOption("objective"),
             optim.method = hydromad.getOption("optim.method"),
             optim.control = hydromad.getOption("optim.control"),
             trace = hydromad.getOption("trace"),
             na.action = na.pass,
             hessian = FALSE)
{
    if (inherits(objective, "formula"))
        objective <- objective[[2]]
    stopifnot(is.language(objective))

    ## TODO: merge this function with hydromad.fit.default?

    optim.control <- modifyList(hydromad.getOption("optim.control"),
                                optim.control)
    ## get data into the right form
    if (is.list(DATA)) DATA <- do.call(ts.intersect, lapply(DATA, as.ts))
    if (!is.ts(DATA)) DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    stopifnot(c("U","Q") %in% colnames(DATA))
    if (is.na(delay))
        delay <- estimateDelay(DATA, plot=FALSE)
    Q <- DATA[,"Q"]
    U <- DATA[,"U"]
    ## initial parameter estimates
    v_s_0 <- 1.0
    ## calculate threshold value of U (switch-over point of volume partitioning)
    ## i.e. value of U where v_s = 0.5
    threshU <- (0.5 / v_s_0) ^ (1 / lambda.init)
    ## split flow record into high/low
    arTime <- autocorrTime(Q, fallTo = 0.02)
    repeat {
        ishigh <- (U >= threshU)
        ishigh[is.na(ishigh)] <- FALSE
        ishigh <- rollmax(ishigh, arTime + delay, na.pad = TRUE, align = "right")
        ishigh[is.na(ishigh)] <- TRUE
        ishigh <- rollmax(ishigh, 2, na.pad = TRUE, align = "left")
        if (sum(ishigh, na.rm = TRUE) < 100) {
            threshU <- mean(U[U > 0], na.rm = TRUE)
            next
        }
        if (sum(!ishigh, na.rm = TRUE) < 100) {
            threshU <- threshU / 2
            next
        }
        break
    }
    lowQ <- Q
    highQ <- Q
    lowQ[ishigh == TRUE] <- NA
    highQ[ishigh == FALSE] <- NA
    ## fit single exponential to low flows
    lowModel <-
        tf.fit(list(Q = lowQ, U = U),
                    order = c(n=1, m=0), delay = delay,
                    method = fit.method,
                    ...)
    if (!inherits(lowModel, "tf"))
        return(lowModel)
    tau_s <- (coef(lowModel, "tau,v"))[["tau_s"]]
    ## fit single exponential to high flows
    highModel <-
        tf.fit(list(Q = highQ, U = U),
                    order = c(n=1, m=0), delay = delay,
                    method = fit.method,
                    ...)
    if (!inherits(highModel, "tf"))
        return(highModel)
    tau_q <- (coef(highModel, "tau,v"))[["tau_s"]]

    ## initial parameter estimates
    init.pars <- c(tau_s = tau_s, tau_q = tau_q, v_s = v_s_0, lambda = lambda.init)
    init.mod <- tf(DATA, init.pars, delay = delay)

    flowpars <- init.pars[which.pars]

    ## set up uh parameter transformations
    uhParTrans <- modifyList(lapply(hydromad.getOption("uhParTrans"), eval),
                             lapply(uhParTrans, eval))
    ## apply transformations
    for (nm in names(flowpars)) {
        ## check for exact match first, then partial match
        fn <- uhParTrans[[nm]]$trans
        if (is.null(fn))
            fn <- uhParTrans[[sub("_.*", "", nm)]]$trans
        if (is.function(fn))
            flowpars[[nm]] <- fn(flowpars[[nm]])
    }
    ## some transformations can induce Inf or -Inf, e.g. log(0)
    make.finite <- function(x)
        if (is.infinite(x)) sign(x) * 100 else x
    flowpars <- lapply(flowpars, make.finite)

    generateModel <- function(pars) {
        flowpp <- as.list(pars)
        ## apply inverse-transform function to each parameter
        for (nm in names(flowpp)) {
            ## check for exact match first, then partial match
            fn <- uhParTrans[[nm]]$inverse
            if (is.null(fn))
                fn <- uhParTrans[[sub("_.*", "", nm)]]$inverse
            if (is.function(fn))
                flowpp[[nm]] <- fn(flowpp[[nm]])
        }
        flowpp <- modifyList(as.list(init.pars), flowpp)
        if (flowpp$tau_q > flowpp$tau_s)
            flowpp[c("tau_q", "tau_s")] <-
                rev(flowpp[c("tau_q", "tau_s")])
        update(init.mod, pars = flowpp)
    }
    xf_optim <- function(pars) {
        #mod <- try(generateModel(pars), silent=TRUE)
        mod <- generateModel(pars)
        if (!isValidModel(mod)) {
            #warning(mod)
            return(NA)
        }
        ## these can be referred to in `objective`
        Q <- observed(mod)
        X <- fitted(mod)
        U <- getU(mod)
        eval(objective)
    }
    ans <- optim(unlist(flowpars), xf_optim, method = optim.method,
                 control = optim.control, hessian = hessian)
    if (hessian) return(ans)
    ## return tf object
    mod <- generateModel(ans$par)
    mod$call <- match.call()
    mod
}




tf.lambda.alt.fit <-
    function(DATA = list(U=, Q=),
             order = c(m = 2, n = 1), ## ignored
             delay = hydromad.getOption("delay"),
             which.pars = c("lambda", "tau_s", "tau_q"), ## tau_q or v_s?
             lambda.init = -0.25,
             fit.method = hydromad.getOption("uh.method"),
             ...,
             uhParTrans = alist(v_s = list(trans = log,
                                inverse = function(x) min(exp(x), 5))),
             objective = hydromad.getOption("objective"),
             optim.method = hydromad.getOption("optim.method"),
             optim.control = hydromad.getOption("optim.control"),
             trace = hydromad.getOption("trace"),
             na.action = na.pass,
             hessian = FALSE)
{
    if (inherits(objective, "formula"))
        objective <- objective[[2]]
    stopifnot(is.language(objective))

    ## TODO: merge this function with hydromad.fit.default?

    optim.control <- modifyList(hydromad.getOption("optim.control"),
                                optim.control)
    ## get data into the right form
    if (is.list(DATA)) DATA <- do.call(ts.intersect, lapply(DATA, as.ts))
    if (!is.ts(DATA)) DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    stopifnot(c("U","Q") %in% colnames(DATA))
    if (is.na(delay))
        delay <- estimateDelay(DATA, plot=FALSE)
    ## split flow record into high/low
    Q <- DATA[,"Q"]
    U <- DATA[,"U"]
    recc <- c(NA, Q) / c(Q, NA)
    recc[recc >= 1] <- NA
    recc[!is.finite(recc)] <- NA
    ## high flow recession rate
    thresh <- mean(Q[Q > quantile(Q, 0.75, na.rm = TRUE)], na.rm = TRUE)
#    highQ <- Q
#    highQ[Q < thresh] <- NA
#    recc <- c(NA, highQ) / c(highQ, NA)
#    recc[recc >= 1] <- NA
    alpha_q <- quantile(recc[Q > thresh], 0.02,
                        na.rm = TRUE, names = FALSE)
    tau_q <- -1 / log(alpha_q)
    ## low flow recession rate
    arTime <- autocorrTime(Q, fallTo = 0.02)
    maskflow <- (is.na(Q) | (Q >= thresh))
    #maskflow[is.na(maskflow)] <- TRUE
    maskflow <- rollmax(maskflow, arTime, na.pad = TRUE, align = "right") #correct?
    maskflow[is.na(maskflow)] <- TRUE
    maskflow <- rollmax(maskflow, 2, na.pad = TRUE, align = "left")
    #minN <- max(hydromad.getOption("warmup") * 1.5, 100)
    #if (sum(maskflow == FALSE, na.rm=TRUE) < minN) {
    #    warning("could not identify low flows (masked by high flow reccessions)")
        alpha_s <- quantile(recc, 0.7, na.rm = TRUE, names = FALSE)
        tau_s <- -1 / log(alpha_s)
    #} else {
    if (FALSE) {
        lowQ <- Q
        lowQ[maskflow == TRUE] <- NA
        ## fit single exponential to low flows
        lowModel <-
            tf.fit(list(Q = lowQ, U = U),
                        order = c(n=1, m=0), delay = delay,
                        method = fit.method,
                        ...)
        if (!inherits(lowModel, "tf"))
            return(lowModel)
        tau_s <- (coef(lowModel, "tau,v"))[["tau_s"]]
        rm(lowModel)
    }
    ## initial parameter estimates
    v_s_0 <- 1.0
    init.pars <- c(tau_s = tau_s, tau_q = tau_q, v_s = v_s_0, lambda = lambda.init)
    init.mod <- tf(DATA, init.pars, delay = delay)

    flowpars <- init.pars[which.pars]

    ## set up uh parameter transformations
    uhParTrans <- modifyList(lapply(hydromad.getOption("uhParTrans"), eval),
                             lapply(uhParTrans, eval))
    ## apply transformations
    for (nm in names(flowpars)) {
        ## check for exact match first, then partial match
        fn <- uhParTrans[[nm]]$trans
        if (is.null(fn))
            fn <- uhParTrans[[sub("_.*", "", nm)]]$trans
        if (is.function(fn))
            flowpars[[nm]] <- fn(flowpars[[nm]])
    }
    ## some transformations can induce Inf or -Inf, e.g. log(0)
    make.finite <- function(x)
        if (is.infinite(x)) sign(x) * 100 else x
    flowpars <- lapply(flowpars, make.finite)

    generateModel <- function(pars) {
        flowpp <- as.list(pars)
        ## apply inverse-transform function to each parameter
        for (nm in names(flowpp)) {
            ## check for exact match first, then partial match
            fn <- uhParTrans[[nm]]$inverse
            if (is.null(fn))
                fn <- uhParTrans[[sub("_.*", "", nm)]]$inverse
            if (is.function(fn))
                flowpp[[nm]] <- fn(flowpp[[nm]])
        }
        flowpp <- modifyList(as.list(init.pars), flowpp)
        if (flowpp$tau_q > flowpp$tau_s)
            flowpp[c("tau_q", "tau_s")] <-
                rev(flowpp[c("tau_q", "tau_s")])
        update(init.mod, pars = flowpp)
    }
    xf_optim <- function(pars) {
        #mod <- try(generateModel(pars), silent=TRUE)
        mod <- generateModel(pars)
        if (!isValidModel(mod)) {
            #warning(mod)
            return(NA)
        }
        ## these can be referred to in `objective`
        Q <- observed(mod)
        X <- fitted(mod)
        U <- getU(mod)
        eval(objective)
    }
    ans <- optim(unlist(flowpars), xf_optim, method = optim.method,
                 control = optim.control, hessian = hessian)
    if (hessian) return(ans)
    ## return tf object
    mod <- generateModel(ans$par)
    mod$call <- match.call()
    mod
}

## TODO: rethink this function

tf.lambda.fit_OLD <-
    function(object,
             ...,
             lambda.init = -0.1,
             which.pars = c("v_s", "lambda", "tau_q"),
             objective = hydromad.getOption("objective"),
             optim.method = hydromad.getOption("optim.method"),
             optim.control = hydromad.getOption("optim.control"),
             hessian = FALSE)
{
    stopifnot(inherits(object, "tf"))
    ## TODO: check that it is a (2, 1) model?
    delay <- object$delay
    ## parse objective function
    if (inherits(objective, "formula")) objective <- objective[[2]]
    stopifnot(is.language(objective))
    theta <- coef(object)#, "tau,v")
    theta <- normalise.tf.coef(theta)
    theta <- as.list(tfParsConvert(theta, "tau,v"))
    ## omit v_q to enforce unit-volume hydrograph
    theta$v_q <- NULL
    if (!("lambda" %in% names(theta)))
        theta$lambda <- lambda.init
    theta <- unlist(theta)

    init.par <- theta[c(which.pars)]
    ## apply logistic transform to lambda
    if ("lambda" %in% which.pars)
        init.par["lambda"] <- qlogis(-init.par["lambda"])

    ih_fn <- function(pars) {
        newpar <- theta
        newpar[names(pars)] <- pars
        ## apply inverse logistic transform to lambda
        if ("lambda" %in% which.pars)
            newpar["lambda"] <- -plogis(newpar["lambda"])
        update(object, pars = newpar)
    }
    ih_optim <- function(pars) {
        model <- try(ih_fn(pars), silent=TRUE)
        if (inherits(model, "try-error")) {
            warning(model)
            return(NA)
        }
        Q <- observed(model)
        X <- fitted(model)
        U <- getU(model)
        eval(objective)
    }
    ans <- optim(init.par, ih_optim, method = optim.method,
                 control = optim.control, hessian = hessian)
    if (hessian) return(ans)
    ## return hydromad object
    ih_fn(ans$par)
}
