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
           hessian = TRUE,
           samples = hydromad.getOption("fit.samples"),
           sampletype = c("latin.hypercube", "random", "all.combinations"),
           multistart = FALSE, initpars = NULL)
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
  \item{hessian}{
    
  }
  \item{samples}{
    number of parameter sets to test.
  }
  \item{sampletype}{
    sampling scheme -- see \code{\link{parameterSets}}.
  }
  \item{multistart}{
%%     ~~Describe \code{multistart} here~~
}
  \item{initpars}{
    initial parameter set.
  }
}
\details{

}
\value{
  the best model from those sampled, according to the given
  \code{objective} function.
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
