\name{defineFeasibleSet}
\alias{defineFeasibleSet}
\alias{defineFeasibleSet.hydromad}
\alias{defineFeasibleSet.default}
\title{
  Compute the set of feasible parameter values meeting some criteria.
}
\description{
  Compute the set of feasible parameter values meeting some
  criteria. These could be used as a feasible set or to estimate
  prediction quantiles according to the 
  GLUE (Generalised Likelihood Uncertainty Estimation) method.
}
\usage{
defineFeasibleSet(x, ...)

\method{defineFeasibleSet}{hydromad}(x, ..., thin = NA)

\method{defineFeasibleSet}{default}(x, objseq, model,
                  frac.within = 0,
                  within.abs = 0.01,
                  within.rel = 0.01,
                  target.coverage = 1,
                  threshold = -Inf,
                  glue.quantiles = NULL,
                  ...)
}
\arguments{
  \item{x}{
    in the \code{hydromad} method this is a \code{hydromad} model object
    which has been run through either \code{\link{fitBySampling}} or
    \code{\link{fitByDream}} to generate a large number of random
    parameter sets with associated objective function scores. 
    In the \code{default} method, \code{x} is a matrix of parameter
    values corresponding to \code{objseq}.
  }
  \item{threshold}{
    value of the objective function (specified in fitting call) used to
    define the feasible set: all parameter sets meeting this criterion
    will be kept. Ignored if \code{NULL}.
  }
  \item{target.coverage}{
    fraction of the observed values to be contained within the ranges of
    simulated values from the feasible set of parameters. 
  }
  \item{eps}{
    observed values must be at least this distance outside the
    uncertainty bounds to count as non-coverage.
  }
  \item{thin}{
    for DREAM results (a Markov Chain Monte Carlo method), the sequences
    should be thinned first to remove autocorrelation and achieve an
    efficient sample. The default is based on the empirical
    autocorrelation but can be specified.
  }
}
\details{
  none yet.
}
\value{
  a modified version of \code{model}, with \code{feasible.set},
  \code{feasible.scores}, \code{feasible.fitted} and \code{coverage}
  values. Can be passed to \code{xyplot}, \code{fitted}, \code{predict},
  \code{update}, \code{print}, etc.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{fitBySampling}},
  \code{\link{fitByDream}}
}
\examples{
data(Queanbeyan)
ts74 <- window(Queanbeyan, start = "1974-01-01", end = "1976-12-01")
mod <- hydromad(ts74, routing = "expuh", rfit = list("inverse", order = c(2,1)))
mod <- update(mod, sma = "cwi",
    tw = c(0,100), f = c(0, 8), loss = c(-0.1, 0.1))
#    routing = "armax", rfit = list("ls", order = c(2,1)))

fit <- fitBySampling(mod, samples = 300, objective = hmadstat("r.sq.sqrt"))

## calculate the set of simulations within 15\% error (or 1 mm/day) 90\% of time
feas <- defineFeasibleSet(fit, frac.within = 0.9, within.rel = 0.15,
                          within.abs = 1)

## view ranges of parameters
feas
## (or)
summary(coef(feas, feasible.set = TRUE))

## plot simulation bounds
xyplot(feas, feasible.bounds = TRUE, cut = 3)

## estimate the 90 percent coverage set, for GLUE.
## calculate 5 percent and 95 percent GLUE quantiles (i.e. weighted).
fitglu <- defineFeasibleSet(fit, target.coverage = 0.9,
  glue.quantiles = c(0.05, 0.95))

## plot simulated GLUE quantiles
xyplot(fitglu, feasible.bounds = TRUE, cut = 3)
xyplot(fitglu, feasible.bounds = TRUE, cut = 3,
       scales = list(y = list(log = TRUE)))
sim <- fitted(fitglu, feasible.bounds = TRUE)
summary(sim$upper - sim$lower)
## coverage of the GLUE quantile simulations
## (not to be confused with the coverage to define feasible set)
mean(sim$lower < observed(fitglu) &
     sim$upper > observed(fitglu))

## simulate on a new data period
newglu <- update(fitglu, newdata = window(Queanbeyan,
                   start = "1980-01-01", end = "1982-01-01"),
                 glue.quantiles = c(0.05, 0.95))
xyplot(newglu, feasible.bounds = TRUE, cut = 3)
sim <- fitted(newglu, feasible.bounds = TRUE)
## coverage of the GLUE quantile simulations
mean(sim$lower < observed(newglu) &
     sim$upper > observed(newglu), na.rm = TRUE)

}
\keyword{ ts }
