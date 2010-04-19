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
    the objective function, as a \code{formula}. See Details.
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
 
  The objective function is given as a \code{\link{formula}} object
  which may refer to the values \code{Q} and \code{X}, representing
  observed and modelled flow, respectively. For more advanced use it may
  also refer to \code{U} (modelled effective rainfall), and more
  generally it may refer to \code{model}, and so may extract other items
  of data, etc.

  The default (unless changed in \code{hydromad.options("objective")} is
  a weighted sum of the R Squared (coefficient of determination) of
  square-root transformed data, and the relative bias.
  
  See \code{\link{summary.hydromad}} for examples of how to specify other
  fit statistics.
  
}
\value{
  the objective function value. It must be numeric and of length one;
  anything else is an error.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \link{hydromad.object},
  \code{\link{hydromad.options}}
}
\examples{

}
\keyword{ utilities }

