\name{objFunVal}
\alias{objFunVal}
\alias{objFunVal.hydromad}
\alias{objFunVal.default}
\title{ Calculate objective function value for a fitted model. }
\description{
  Calculate objective function value for a fitted model.
}
\usage{
objFunVal(x, objective, ...)

\method{objFunVal}{hydromad}(x, objective = hydromad.getOption("objective"),
          ..., all = FALSE, nan.ok = FALSE)

\method{objFunVal}{default}(x, objective = hydromad.getOption("objective"),
          ..., nan.ok = FALSE)
}
\arguments{
  \item{x}{
    object from which to calculate stats. For the \code{hydromad} method,
    this should be a fitted \code{hydromad} model object, i.e. it must
    have specific parameter values, not ranges. For the default method,
    this should be a matrix-like object with named columns \code{Q},
    \code{X} and optionally \code{U} and \code{P}.
  }
  \item{objective}{
    the objective function, or a list of objective functions. See Details.
  }
  \item{...}{
    ignored.
  }
  \item{all}{
    passed to \code{fitted} and \code{observed}.
  }
  \item{nan.ok}{
    by default, an error is thrown if the result is \code{NaN}. Set this
    argument to \code{TRUE} to suppress the error.
  }
}
\details{
 
  The objective function is given as a function with arguments \code{Q},
  \code{X} and \code{\dots}, and optionally other arguments.  \code{Q}
  and \code{X} represent observed and modelled flow, respectively.
  It should return a single numeric value.
  
  For more advanced use it may also refer to arguments \code{U}
  (modelled effective rainfall), or \code{P} (observed rainfall), and
  more generally it may refer to \code{model}, and so may extract other
  items of data, parameters, etc.

  The default (unless changed in \code{hydromad.options("objective")}) is
  a weighted sum of the R Squared (coefficient of determination) of
  square-root transformed data, and the relative bias.
  
  See \code{\link{hydromad.stats}} for examples of how to specify other
  fit statistics.
}
\value{
  the objective function value, or a list of objective function
  values. They must be numeric and of length one; anything else is an
  error.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad.stats}},
  \link{hydromad.object}
}
\examples{
dat <- data.frame(Q = rnorm(10), X = rnorm(10))
objFunVal(dat, hmadstat("RMSE"))
}
\keyword{ utilities }

