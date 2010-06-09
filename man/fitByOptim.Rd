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
           method = hydromad.getOption("optim.method"), control = list(),
           samples = hydromad.getOption("fit.samples"),
           sampletype = c("latin.hypercube", "random", "all.combinations"),
           initpars = NULL, multistart = FALSE,
           vcov = FALSE, hessian = vcov)

fitByOptim1(MODEL, objective = hydromad.getOption("objective"),
            tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{MODEL}{
    a model specification created by \code{\link{hydromad}}. It should
    not be fully specified, i.e one or more parameters should be defined
    by \emph{ranges} of values rather than exact values.
  }
  \item{objective}{
    objective function to minimise, given as a \code{function(Q, X, ...)}.
    See \code{\link{objFunVal}}.
  }
  \item{method, control}{
    optimisation algorithm and settings. See \code{\link{optim}}. The
    default is \code{"Nelder-Mead"} (a simplex algorithm).
    An additional method \code{"PORT"} can be given, which calls
    \code{\link{nlminb}}. 
  }
  \item{samples}{
    number of parameter sets to test.
  }
  \item{sampletype}{
    sampling scheme -- see \code{\link{parameterSets}}.
  }
  \item{initpars}{
    initial parameter set.
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
  \item{tol}{
    convergence tolerance, see \code{\link{nlminb}}. 
  }
}
\details{
  See \code{\link{optim}} for a brief description of the available
  algorithms, including references to literature.
}
\value{
  the best model from those sampled, according to the given
  \code{objective} function. Also, these extra elements are inserted:
  \item{fit.result}{
    the result from \code{\link{optim}}.
  }
  \item{objective}{
    the \code{objective} function used.
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
  \code{\link{fitBySampling}},
  \code{\link{optim}},
  \code{\link{objFunVal}}
}
\examples{
data(Cotter)
x <- Cotter[1:1000]

## IHACRES CWI model with power law unit hydrograph
modx <- hydromad(x, sma = "cwi", routing = "powuh")
modx

foo <- fitByOptim(modx, vcov = TRUE)

summary(foo)

## parameter correlation matrix (dodgy!) with symbols
symnum(cov2cor(vcov(foo)))

## return value from optim (with 'objseq' added):
str(foo$fit.result)

## plot objective function value convergence over time
ppp <- foo$fit.result$objseq
convergenceTs <- ts(ppp, end = foo$funevals,
  deltat = foo$funevals / length(ppp))
xyplot(convergenceTs, type = "b",
  xlab = "function evaluations", ylab = "objective fn. value")
}
\keyword{ optimization }
