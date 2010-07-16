

glue.hydromad <-
    function(model, ...,
             target.coverage = 0.9, eps = 0.01,
             max.n = 250, thin = NA)
{
    stopifnot(inherits(model$fit.result, "dream"))
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
    mcmc1 <- as.matrix(mcmc)
    logp <- as.vector(logp)
    ordlik <- order(logp, decreasing = TRUE)
    ## extract observed data to check coverage of uncertainty bounds
    obs <- observed(model)
    ## best model as starting point
    xsim <- fitted(update(model, newpars = mcmc1[ordlik[1],]))
    ## uncertainty bounds
    sim.lower <- sim.upper <- xsim
    cover <- 0
    max.n <- min(max.n, nrow(mcmc1))
    for (i in 2:max.n) {
        xsim <- fitted(update(model, newpars = mcmc1[ordlik[i],]))
        ## uncertainty bounds
        sim.lower <- pmin(sim.lower, xsim)
        sim.upper <- pmax(sim.upper, xsim)
        ## check coverage interval
        cover <- mean((sim.lower - eps < obs) & (obs < sim.upper + eps))
        if (cover >= target.coverage) break
    }
    list(sim = cbind(obs = obs, lower = sim.lower, upper = sim.upper),
         coverage = cover, thin = thin, 
         pars = mcmc1[ordlik[1:i],])
}
