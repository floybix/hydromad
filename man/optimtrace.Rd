\name{optimtrace}
\alias{optimtrace}
\alias{optimtrace.default}
\title{ Extract objective function value series from optimisation results }
\description{
  Extract objective function value series from optimisation results. It
  is assumed that the results are from a minimisation problem.
}
\usage{
optimtrace(object, ...)
}
\arguments{
  \item{object}{
    fitted model object. Currently only \code{\link{hydromad}} objects
    are supported.
  }
  \item{\dots}{
    further arguments which may be supported for particular optimisation
    algorithms. These are:
    \describe{
      \item{\code{raw}}{
        if \code{TRUE}, then raw function values are returned, which may
        consist of multiple series in parallel, and may jump up and
        down. Otherwise (the default), the the best value so far at each
        time step is returned.
      }
      \item{\code{objective}}{
        a function or formula to calculate from \code{object} to
        construct the objective function value series.
        It is passed to \code{\link{objFunVal}}.
        Currently this is only supported for \code{dream} results.
      }
    }
  }
}
\value{
  a \code{zoo} series where the time index represents number of function
  evaluations. 
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{fitBySCE}},
  \code{\link{fitByDE}},
  \code{\link{fitByDream}},
  \code{\link{fitByOptim}}
}
\examples{

}
\keyword{ optimization }
