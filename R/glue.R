

glue.hydromad <-
    function(model,
             threshold.value = NULL,
             target.coverage = 0.9, eps = 0.01,
             thin = NA)
{
    if (inherits(model$fit.result, "dream")) {
        result <- model$fit.result
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
        result <- model$fit.result
        objseq <- result$objseq
        psets <- result$psets
        if (is.null(objseq) || is.null(psets)){
            stop(paste("this function only works on results from fitByDream() or fitBySampling()",
                       "(or other methods generating fit.result$objseq and fit.result$psets)"))
        }
    }
    glue.hydromad.default(model, objseq = objseq, psets = psets,
                          threshold = threshold.value,
                          target.coverage = target.coverage, eps = eps)
}


glue.hydromad.default <-
    function(model, objseq, psets,
             threshold = NULL,
             target.coverage = NULL, eps = 0.01)
{
    ## order parameter sets by objective function value
    ordlik <- order(objseq, decreasing = TRUE)
    objseq <- objseq[ordlik]
    psets <- psets[ordlik,,drop=FALSE]
    ## extract observed data to check coverage of uncertainty bounds
    obs <- observed(model)
    ## best model as starting point
    sim.lower <- sim.upper <- fitted(model, all = TRUE)
    cover <- 0
    for (i in 2:length(objseq)) {
        ## check coverage interval
        isinside <- (sim.lower - eps < obs) & (obs < sim.upper + eps)
        cover <- mean(isinside[-(1:model$warmup)], na.rm = TRUE)
        if (!is.null(target.coverage) && (cover >= target.coverage))
            break
        thisPars <- as.list(psets[i,])
        thisVal <- objseq[i]
        if (!is.finite(thisVal))
            break
        ## don't expand feasible set beyond the threshold (if given)
        if (!is.null(threshold) && (thisVal < threshold))
            break
        xsim <- fitted(update(model, newpars = thisPars), all = TRUE)
        sim.lower <- pmin(sim.lower, xsim)
        sim.upper <- pmax(sim.upper, xsim)
    }
    ok <- 1:(i-1)
    model$feasible.set <- as.matrix(psets[ok,,drop=FALSE])
    model$feasible.scores <- objseq[ok]
    model$feasible.fitted <- cbind(lower = sim.lower,
                                   upper = sim.upper)
    model$feasible.coverage <- cover
    model
}
