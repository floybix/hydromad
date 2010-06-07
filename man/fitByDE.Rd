\name{fitByDE}
\alias{fitByDE}
\title{
  Fit a hydromad model using the DE (Differential Evolution) algorithm.
}
\description{
  Fit a hydromad model using the DE (Differential Evolution) algorithm.
}
\usage{
fitByDE(MODEL, 
        objective = hydromad.getOption("objective"),
        control = hydromad.getOption("de.control"))
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
  \item{control}{
    settings for the DE algorithm. See \code{\link[DEoptim]{DEoptim.control}}.
  }
}
\value{
  the best model from those sampled, according to the given
  \code{objective} function. Also, these extra elements are inserted:
  \item{fit.result}{
    the result from \code{\link[DEoptim]{DEoptim}}.
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
  \code{\link[DEoptim]{DEoptim}},
  \code{\link{objFunVal}}
}
\examples{
library("DEoptim")

data(Cotter)
x <- Cotter[1:1000]

## IHACRES CWI model with power law unit hydrograph
modx <- hydromad(x, sma = "cwi", routing = "powuh")
modx

foo <- fitByDE(modx, control = DEoptim.control(itermax = 10))

summary(foo)

## return value from DE:
str(foo$fit.result)

## plot objective function value convergence over time
ppp <- foo$fit.result$member$bestvalit
convergenceTs <- ts(foo$fit.result$member$bestvalit,
  end = foo$funevals, deltat = foo$funevals / length(ppp))
xyplot(convergenceTs, type = "b",
  xlab = "function evaluations", ylab = "objective fn. value")
}
\keyword{ optimization }
