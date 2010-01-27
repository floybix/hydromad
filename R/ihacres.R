## ihacreslab: rainfall-runoff hydrology models and tools
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


### EXAMPLE SESSION ###
#
# ts70s <- SalmonBrook['1970::1980']
# foo <- ihacres(ts70s, sma = "cwi",
#                tw = range(2, 128),
#                f = range(0.1, 3),
#                routing = "expuh", warmup = 50,
#                rfit = list("ls", order = c(1,0), delay = 0))

# foo <- ihacres(ts70s, sma = "cwi",
#                routing = "expuh", tau_s = 30, delay = 0)

#
# ### NOTE absorbScale()
#
# fitBySampling(foo, samples = 64)
#
# fitByOptim(foo, samples = 50)
#
# foo <- ihacres(ts70s, sma = "cwi", l = range(0,200), tau_s = 10)
#
# foo <- fitBySCE(foo)
#
# foo <- update(foo, rfit = list("sriv", order = c(1,1)))
#
# mcmcByDRAM(foo)
# mcmcByDream(foo)


# ## define composite model with multiple rain gauge inputs
# ## (model parameters are common across all)
# ## (weights for subcatchments become [n-1] extra parameters)

# ## Question: should U from each subcatchment be summed before routing;
# ##           or should U from each be routed, and each Q summed?
# ## -- Probably the latter, so that lambda and loss make sense.
# ##    However, this means we will have to handle all the tf() methods:
# ##    fitted, predict, summary, etc
# ##    Unless - generic tf can handle multiple U inputs and sum the results...
# ##             but will need weights for summing.
# ##             (only works when tf parameters are same across all)
# ##    Alternative: do routing all within SMA; null tf() only.

####### TODO: if everything is passed through to tf(pars=),
#######       can it detect length(pars) == 0 ???

# moo <- multip.ihacres(mults[,1], PDATA = mults[,-(1)],
#                       sma = "cwi", f = 0, rfit = list("ls"))
# fitByOptim(moo, samples = 50)
#
# ## specified weights (area-weighted) -- constrain sum(weights) == 1
# moo2 <- update(moo, weight.1 = 0.3, weight.2 = 0.5)
# fitByOptim(moo2)


multi.cmd.sim <-
    function(DATA = zoo(),
             d, f, e, M_0 = d/2,
             weight)
{}

multi.sma.sim.fun <-
    function(sma = "cwi")
{}

multi.cwi.sim <-
    function(DATA = zoo(),
             tw, f = 0, c,
             l = 0, p = 1,
             t_ref = ihacres.getOption("cwi")$t_ref,
             s_0 = 0,
             ...)
{
    DATA <- as.zoo(DATA)
    whichP <- ((colnames(DATA) %in% c("Q","E")) == FALSE)
    PDATA <- DATA[,whichP]
    DATA <- DATA[,!whichP]
    #dots <- list(...)
    #whichW <- grepl("^weight", names(dots))
    #weights <- dots[whichW]
    #sma.pars <- dots[!whichW]

    ## how to pass only relevant arguments to this (SMA) function?
    ## if everything goes through dots, what about UH parameters?

    weights <- list(...)
    weights <- weights[grep("^w", names(weights))]
    sma.pars <- list(tw = tw, f = f, c = c, l = l, p = p, t_ref = t_ref, s_0 = s_0)
    U <- PDATA[,1] * 0
    uList <- list()
    for (i in 1:NCOL(PDATA)) {
        uCall <- quote(cwi.sim(merge(DATA, P = PDATA[,i])))
        uCall <- as.call(c(as.list(uCall),
                           sma.pars))
        iU <- eval(uCall)
        U <- U + iU * weights[[i]]
    }
    U
}

absorbScale.multi.cwi <- absorbScale.cwi

multip.ihacres <-
    function(DATA = zoo(),
             PDATA = zoo(),
             ...,
             sma = "cwi",
             routing = "expuh")
{
    DATA <- as.zoo(DATA)
    PDATA <- as.zoo(PDATA)
    nRain <- NCOL(PDATA)


    multi.sma.fun <- function(DATA, ..., weights)
    {

    }
    mod0 <- ihacres(DATA, ..., class = class,
                    uh = FALSE)


}

#
#
# ihacres.options(sacramento = list(lzpr = range(0, 100), ...))
# bar <- ihacres(ts70s, sma = "sacramento",
#                routing = "expuh", rfit = "inverse")
#
# bar <- fitBySCE(bar)
#
# mcmcByDream(bar)
#
#
# bucket.sim <- function(DATA, capacity, evap, w=0) {
#   P <- DATA$P
#   U <- P
#   for (t in 1:NROW(P)) {
#     w <- w + P[t]
#     U[t] <- max(0, w - capacity)
#     w <- max(0, w - U[t] - evap)
#   }
#   return(U)
# }
# baz <- ihacres(ts70s, sma = "bucket",
#          capacity = range(1, 100), evap = range(0, 10), rfit = "ls")
# fitByOptim(baz, samples = 50, method = "Nelder-Mead")
# mcmcByDRAM(baz)
#
#
########################


ihacres <-
    function(DATA = zoo(),
             ...,
             sma = ihacres.getOption("sma"),
             routing = ihacres.getOption("routing"),
             rfit = NULL,
             weights = NULL,
             warmup = ihacres.getOption("warmup"))
{
    ## TODO: keep zoo attributes?
    DATA <- as.ts(DATA)
    ## dots `...` may contain arguments for sma and/or routing.
    ## parlist stores these -- and may be ranges (!isFullySpecified).
    ## take defaults from ihacres.options()
    parlist <- list()
    if (is.character(routing)) {
        parlist <-
            modifyList(parlist, as.list(ihacres.getOption(routing)))
    }
    if (is.character(sma)) {
        parlist <-
            modifyList(parlist, as.list(ihacres.getOption(sma)))
    }
    parlist <- modifyList(parlist, list(...))
    ## create the model object
    obj <- list(call = match.call())
    class(obj) <- unique(c(sma, "ihacres"))
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

isFullySpecified <-
    function(object, which = c("both", "sma", "routing"))
{
    !is.list(coef(object, which = which, warn = FALSE))
}

coef.ihacres <-
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
    if ("..." %in% sma.argnames) {
        forSMA <- forSMA | unmatched
    }
    if ("..." %in% r.argnames) {
        forRouting <- forRouting | unmatched
    }
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

fitted.ihacres <- function(object, ..., all = FALSE)
{
    tmp <- object$fitted.values
    if (is.null(object$routing))
        tmp <- object$U
    if (all) return(tmp)
    stripWarmup(tmp, object$warmup)
}

residuals.ihacres <- function(object, ..., all = FALSE)
{
    tmp <- (object$data[,"Q"] - fitted(object, all = TRUE))
    if (all) return(tmp)
    stripWarmup(tmp, object$warmup)
}

observed.ihacres <- function(object, ..., all = FALSE)
{
    ## observed.default will work, but this may be slightly faster
    tmp <- object$data[,"Q"]
    if (all) return(tmp)
    stripWarmup(tmp, object$warmup)
}

vcov.ihacres <- function(object, ...)
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

deviance.ihacres <- stats:::deviance.lm

getU <- function(object, ...)
    UseMethod("getU")

getU.ihacres <- function(object, ..., all = FALSE)
{
    tmp <- object$U
    if (all) return(tmp)
    stripWarmup(tmp, object$warmup)
}

getP <- function(object, ...)
    UseMethod("getP")

getP.ihacres <- function(object, ..., all = FALSE)
{
    if (!("P" %in% colnames(object$data)))
        return(NULL)
    tmp <- object$data[,"P"]
    if (all) return(tmp)
    stripWarmup(tmp, object$warmup)
}

print.ihacres <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    rx <- x$data
    routingStr <- if (is.null(x$routing)) "NULL" else toString(x$routing)
    cat(paste("\nIHACRES model of (SMA) class \"", class(x)[1],
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
                print(t(sapply(parlist,
                               function(p) c(lower = min(p), upper = max(p)))),
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

isValidModel.ihacres <- function(object, ...)
    is.numeric(fitted(object, all = TRUE))
