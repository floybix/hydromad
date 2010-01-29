\name{tryModelOrders}
\alias{tryModelOrders}
\title{
  Compare calibrations with different order models.
}
\description{
  Compare calibrations with different order models.
}
\usage{
tryModelOrders(expr, n = 0:3, m = 0:2,
               delay = ihacres.getOption("delay"),
               verbose = hydromad.getOption("trace"))
}
\arguments{
  \item{expr}{
    an expression to calibrate an hydromad type model.
  }
  \item{n}{
    the set of values of n to try (order of the auto-regressive component).
  }
  \item{m}{
    the set of values of m to try (order of the moving-average component).
  }
  \item{delay}{
    the set of delays to try.
  }
  \item{verbose}{
    show detailed tracing output.
  }
}
\details{

}
\value{
  a list of model objects, of class \code{hydromad.\link{runlist}}.
}
\references{
  P C Young?
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{tf}},
  \code{\link{tf.sriv.fit}}
}
\examples{

}
\keyword{ optimize }
