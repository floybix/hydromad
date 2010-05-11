\name{fitByDream}
\alias{fitByDream}
\title{
  Fit a hydromad model using the DREAM (DiffeRential Evolution Adaptive Metropolis) algorithm.
}
\description{
  Fit a hydromad model using the DE (DiffeRential Evolution Adaptive Metropolis) algorithm.
}
\usage{
fitByDream(MODEL, loglik = ~ -0.5 * sum((Q-X)^2),
           control = hydromad.getOption("dream.control"),
           vcov = TRUE)
}
\arguments{
  \item{MODEL}{
    a model specification created by \code{\link{hydromad}}. It should
    not be fully specified, i.e one or more parameters should be defined
    by \emph{ranges} of values rather than exact values.
  }
  \item{loglik}{
    log-likelihood function (log of the posterior probability density),
    given as a \code{formula} involving 
    \code{Q}, \code{X} and/or \code{U}. See \code{\link{objFunVal}}.
  }
  \item{control}{
    settings for the DREAM algorithm. See \code{\link[dream]{dream}}.
  }
  \item{vcov}{
    if \code{vcov = TRUE}, the parameter variance-covariance matrix
    will be estimated from the last half of the sequences.
    It can be extract using \code{\link{vcov}}.
  }
}
\value{
  the best model from those sampled, according to the given
  \code{loglik} function. Also, these extra elements are inserted:
  \item{fit.result}{
    the result from \code{\link[dream]{dream}}.
  }
  \item{objective}{
    the \code{loglik} function used.
  }
  \item{funevals}{
    total number of evaluations of the model simulation function.
  }
  \item{timing}{
    timing vector as returned by \code{system.time}.
  }
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link[dream]{dream}},
  \code{\link{objFunVal}}
}
\examples{
if (require("dream")) {

  data(Cotter)
  x <- Cotter[1:1000]

  ## IHACRES CWI model with power law unit hydrograph
  modx <- hydromad(x, sma = "cwi", routing = "powuh")
  modx
  
  foo <- fitByDream(modx, control = list(maxtime = 30))

  summary(foo)

  ## parameter correlation matrix with symbols
  symnum(cov2cor(vcov(foo)))

  ## return value from dream:
  str(foo$fit.result)

  ppp <- foo$fit.result$hist.logp
  convergenceTs <- ts(-ppp, deltat = foo$funevals / nrow(ppp))
  xyplot(window(convergenceTs, start = 50),
    superpose = TRUE, auto.key = FALSE,
    xlab = "function evaluations", ylab = "neg. log likelihood")

  ## MCMC diagnostics:
  summary(foo$fit.result)
  plot(foo$fit.result)

}
}
\keyword{ optimization }
