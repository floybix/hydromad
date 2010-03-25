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

\S3method{summary}{runlist}(object, ..., FUN = summary)
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
    a \code{runlist}: a list of fitted model objects.
  }
  \item{FUN}{
    function returning one or more named numeric values. Any returned
    values other than single numeric values will be ignored.
  }
}
\details{
  Note that the \code{coef} method just calls \code{summary(..., FUN = coef)}.
}
\value{
  \code{runlist} and \code{as.runlist} return a list of class \code{runlist}
  and also \var{class}\code{.runlist}, where \var{class} is the first
  class of the first element of the list.

  \code{summary} and \code{coef} return a data frame, with rows
  for each element of \code{object} and columns for each named item
  returned by \code{FUN}. Any missing items will filled in with \code{NA}.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{xyplot.runlist}},
  \code{\link{coef.hydromad}},
  \code{\link{summary.hydromad}}
}
\examples{

}
\keyword{ utilities }


