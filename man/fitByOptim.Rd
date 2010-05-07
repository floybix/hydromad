\name{fitByOptim}
\alias{fitByOptim}
\title{
  Fit a hydromad model using a general-purpose optimiser.
}
\description{
  Fit a hydromad model using a general-purpose optimiser.
}
\usage{
fitByOptim(MODEL, objective = hydromad.getOption("objective"),
           method = hydromad.getOption("optim.method"),
           control = hydromad.getOption("optim.control"),
           samples = hydromad.getOption("fit.samples"),
           sampletype = c("latin.hypercube", "random", "all.combinations"),
           multistart = FALSE,
           vcov = FALSE, hessian = vcov, initpars = NULL)
}
\arguments{
  \item{MODEL}{
    a model specification created by \code{\link{hydromad}}. It should
    not be fully specified, i.e one or more parameters should be defined
    by \emph{ranges} of values rather than exact values.
  }
  \item{objective}{
    objective function to minimise, given as a \code{formula} involving
    \code{Q}, \code{X} and/or \code{U}. See \code{\link{objFunVal}}.
  }
  \item{method, control}{
    optimisation algorithm and settings. See \code{\link{optim}}. The
    default, "Nelder-Mead" (a simplex algorithm), is recommended.
  }
  \item{samples}{
    number of parameter sets to test.
  }
  \item{sampletype}{
    sampling scheme -- see \code{\link{parameterSets}}.
  }
  \item{multistart}{
    if this is \code{TRUE}, then each of the initial parameter samples
    (based on \code{samples} and \code{sampletype}) are used as a
    starting value for a separate optimisation run, and the best result
    is kept.
  }
  \item{vcov, hessian}{
    if \code{vcov = TRUE}, the parameter variance-covariance matrix
    will be estimated as the inverse of
    the hessian matrix returned by \code{optim()}. This may be a very
    poor estimate! It can be extract using \code{\link{vcov}}.
  }
  \item{initpars}{
    initial parameter set.
  }
}
\details{

}
\value{
  the best model from those sampled, according to the given
  \code{objective} function. The result from \code{\link{optim}} is
  stored as item \code{fit.result}.
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{fitBySampling}},
  \code{\link{optim}},
  \code{\link{objFunVal}}
}
\examples{

}
\keyword{ optimization }
