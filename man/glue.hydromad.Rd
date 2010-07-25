\name{glue.hydromad}
\alias{glue.hydromad}
\alias{glue.hydromad.default}
\title{
  Extract a feasible parameter set, like GLUE method.
}
\description{
  Extract a feasible parameter set, like GLUE method.
}
\usage{
glue.hydromad(model,
              threshold.value = NULL,
              target.coverage = 0.9, eps = 0.01,
              thin = NA)
}
\arguments{
  \item{model}{
    a \code{hydromad} object which has been run through either
    \code{\link{fitBySampling}} or \code{\link{fitByDream}} to generate
    a large number of random parameter sets with associated objective
    function scores. 
  }
  \item{threshold.value}{
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
mod <- hydromad(ts74, sma = "cwi",
    routing = "armax", rfit = list("ls", order = c(2,1)))
fit <- fitBySampling(mod, samples = 100)
fitglu <- glue.hydromad(fit, target.coverage = 0.9)
fitglu
coef(fitglu, feasible.set = TRUE)
xyplot(fitglu, feasible.bounds = TRUE)
}
\keyword{ ts }
