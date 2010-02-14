\name{scalar}
\alias{scalar}
\alias{scalar.sim}
\title{ Simple constant runoff proportion }
\description{
  Simple constant runoff proportion.
}
\usage{
scalar.sim(DATA, c, return_state = FALSE)
}
\arguments{
  \item{DATA}{time-series-like object with a column P (precipitation).}
  \item{c}{fraction of rainfall that becomes effective.}
  \item{return_state}{ignored.}
}
\value{
  the simulated effective rainfall, a time series of the same length as
  the input series.  
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}(sma = "scalar")} to work with models as objects (recommended).
}
\examples{
## view default parameter ranges:
str(hydromad.options("scalar"))

data(Cotter)
sca0 <- hydromad(Cotter[1:500], sma = "scalar", routing = "expuh")
sca0
}
\keyword{models}
