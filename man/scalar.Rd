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

data(HydroTestData)
mod0 <- hydromad(HydroTestData, sma = "scalar", routing = "expuh")
mod0

## simulate with some arbitrary parameter values
testQ <- predict(update(mod0, c = 0.5, tau_s = 10))
xyplot(ts.union(HydroTestData[,1:2], scalar.Q = testQ))
}
\keyword{models}
