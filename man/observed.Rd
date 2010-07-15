\name{observed}
\alias{observed}
\alias{observed.default}
\title{ Observed data values }
\description{
	Like \code{\link{fitted}}, but return the original "observed" data
	rather than the modelled data.
}
\usage{
observed(object, ...)

\method{observed}{default}(object, ...)
}
\arguments{
  \item{object}{ a model object for which to return the observed values. }
  \item{\dots}{ other arguments. }
}
\details{
	The default method should work for any object with methods for \code{residuals} and \code{fitted}.
}
\value{
  the result of \code{fitted(object) + residuals(object)}.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{ \code{\link{residuals}} }
\examples{
x <- 1:10
mymodel <- list(fitted.values = x + rnorm(10))
mymodel$residuals <- x - fitted(mymodel)
observed(mymodel)
}
\keyword{ ts }

