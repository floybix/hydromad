

defineFeasibleSet <- function(x, ...)
    UseMethod("defineFeasibleSet")

defineFeasibleSet.hydromad <-
    function(x, ..., thin = NA)
{
    x <- update(x, feasible.set = NULL)
    if (inherits(x$fit.result, "dream")) {
        result <- x$fit.result
        ## extract last half of Sequences (assumed to have converged)
        mcmc <- window(result)
        logp <- as.mcmc(window(result$hist.logp, start = start(mcmc)))
        ## thin the sequences to remove autocorrelation (more efficient samples)
        if (is.na(thin)) {
            effsz <- effectiveSize(mcmc)
            thin <- 2 * ceiling(max(nrow(mcmc[[1]]) / effsz))
        }
        mcmc <- window(mcmc, thin = thin)
        logp <- window(logp, thin = thin)
        ## stack all parallel MCMC chains into one
        psets <- as.matrix(mcmc)
        objseq <- as.vector(logp)
    } else {
        result <- x$fit.result
        objseq <- result$objseq
        psets <- result$psets
        if (is.null(objseq) || is.null(psets)){
            stop(paste("this function only works on results from fitByDream() or fitBySampling()",
                       "(or other methods generating fit.result$objseq and fit.result$psets)"))
        }
    }
    defineFeasibleSet.default(psets, objseq = objseq, model = x, ...)
}

defineFeasibleSet.default <-
    function(x, objseq, model,
             frac.within = 0,
             within.abs = 0.01,
             within.rel = 0.01,
             target.coverage = 1,
             threshold = -Inf,
             glue.quantiles = NULL,
             ...)
{
    psets <- x
    model <- update(model, feasible.set = NULL)
    ## order parameter sets by objective function value
    ordlik <- order(objseq, decreasing = TRUE)
    objseq <- objseq[ordlik]
    psets <- psets[ordlik,,drop=FALSE]
    ## extract observed data to check coverage of uncertainty bounds
    obs <- observed(model, all = TRUE)
    ## tolerances for each time step
    obs.tol <- pmax(abs(obs * within.rel), within.abs)

    sims <- NULL
    sim.lower <- sim.upper <- NULL
    ok <- rep(FALSE, length(obs))
    
    for (i in 2:length(objseq)) {
        thisPars <- as.list(psets[i,])
        thisVal <- objseq[i]
        if (!is.finite(thisVal))
            break
        ## don't expand feasible set beyond the threshold objective (if given)
        if (thisVal < threshold)
            break
        ok[i] <- TRUE
        ## simulate
        xsim <- fitted(update(model, newpars = thisPars), all = TRUE)
        ## check what fraction of this simulation is within given tolerances
        if (frac.within > 0) {
            this.ok <- abs(xsim - obs < obs.tol)
            this.frac.within <- mean(this.ok[-(1:model$warmup)], na.rm = TRUE)
            if (this.frac.within < frac.within) {
                ok[i] <- FALSE
                next
            }
        }
        if (is.null(sim.lower)) {
            sim.lower <- sim.upper <- xsim
        }
        sim.lower <- pmin(sim.lower, xsim)
        sim.upper <- pmax(sim.upper, xsim)
        if (!is.null(glue.quantiles)) {
            sims <- cbind(sims, coredata(xsim))
        }
        ## check coverage of cumulative set of simulations
        if (target.coverage < 1) {
            isinside <- (sim.lower - within.abs < obs) & (obs < sim.upper + within.abs)
            cover <- mean(isinside[-(1:model$warmup)], na.rm = TRUE)
            if (cover > target.coverage)
                break
        }
    }
    model$feasible.set <- as.matrix(psets[ok,,drop=FALSE])
    model$feasible.scores <- objseq[ok]
    if (!is.null(glue.quantiles)) {
        glue.quantiles <- range(glue.quantiles)
        weights <- model$feasible.scores - min(model$feasible.scores, na.rm = TRUE)
        weights <- weights / sum(weights, na.rm = TRUE)
        bounds <-
            t(apply(sims, 1, safe.wtd.quantile, weights = weights,
                    probs = glue.quantiles, normwt = TRUE))
        colnames(bounds) <- c("lower", "upper")
        model$feasible.fitted <- zoo(bounds, time(obs))
    } else {
        model$feasible.fitted <- cbind(lower = sim.lower,
                                       upper = sim.upper)
    }
    model
}

safe.wtd.quantile <- function(x, ..., probs) {
    if (any(is.finite(x)))
        wtd.quantile(x, ..., probs = probs)
    else rep(NA_real_, length(probs))
}
