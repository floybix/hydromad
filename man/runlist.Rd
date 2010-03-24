\name{runlist}
\alias{runlist}
\alias{as.runlist}
\alias{[.runlist}
\alias{coef.runlist}
\alias{summary.runlist}
\alias{print.summary.runlist}
\title{ Work with a set of model runs }
\description{
  A \code{runlist} object is simply a list of model objects.
}
\usage{
runlist(...)

as.runlist(x, ...)

\S3method{coef}{runlist}(object, ...)

\S3method{summary}{runlist}(object, ...)
}
\arguments{
  \item{\dots}{
    for \code{runlist}, a named list of model objects,
    specified directly as in \code{\link{list}}.

    In other cases, arguments are passed on to the generic functions. 
  }
  \item{x}{
    a simple \code{list}.
  }
  \item{object}{
    a \code{runlist}.
  }
}
\details{

}
\value{

}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{xyplot.runlist}}
}
\examples{

}
\keyword{ utilities }


