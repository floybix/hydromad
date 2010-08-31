

defineFeasibleSet <- function(x, ...)
    UseMethod("defineFeasibleSet")

defineFeasibleSet.hydromad <-
    function(x, ..., thin = NA)
{
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

## TODO

defineFeasibleSet.default <-
    function(x, objseq, model,
             within.abs = 0.01,
             within.rel = 0.01,
             frac.within = 0,
             target.coverage = 1,
             threshold = -Inf,
             glue.quantiles = NULL)
{
    psets <- x
    ## order parameter sets by objective function value
    ordlik <- order(objseq, decreasing = TRUE)
    objseq <- objseq[ordlik]
    psets <- psets[ordlik,,drop=FALSE]
    ## extract observed data to check coverage of uncertainty bounds
    obs <- observed(model, all = TRUE)
    ## best model as starting point
    bestmodel <- update(model, newpars = as.list(psets[1,]))
    sim.lower <- sim.upper <- xsim <- fitted(bestmodel, all = TRUE)
    cover <- 0
    for (i in 2:length(objseq)) {
        thisPars <- as.list(psets[i,])
        thisVal <- objseq[i]
        if (!is.finite(thisVal))
            break
        ## don't expand feasible set beyond the threshold (if given)
        if (thisVal < threshold)
            break
        xsim <- fitted(update(model, newpars = thisPars), all = TRUE)
        sim.lower.tmp <- pmin(sim.lower, xsim)
        sim.upper.tmp <- pmax(sim.upper, xsim)
        ## check what fraction of this simulation is within given tolerances
        this.eps <- pmax(abs(xsim * within.rel), within.abs)
        this.ok <- (xsim - this.eps < obs) & (obs < xsim + this.eps)
        this.within <- mean(this.ok[-(1:model$warmup)], na.rm = TRUE)
        if (this.within < frac.within)
            break
        ## check coverage of cumulative set of simulations
        isinside <- (sim.lower.tmp - within.abs < obs) & (obs < sim.upper.tmp + within.abs)
        cover <- mean(isinside[-(1:model$warmup)], na.rm = TRUE)
        if (cover > target.coverage)
            break
        sim.lower <- sim.lower.tmp
        sim.upper <- sim.upper.tmp
    }
    ok <- 1:(i-1)
    model$feasible.set <- as.matrix(psets[ok,,drop=FALSE])
    model$feasible.scores <- objseq[ok]
    model$feasible.fitted <- cbind(lower = sim.lower,
                                   upper = sim.upper)
    model
}
