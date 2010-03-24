\name{dbm}
\alias{dbm}
\alias{dbm.sim}
\title{ Typical initial model used in Data-Based Mechanistic modelling. }
\description{
  Typical initial model used in Data-Based Mechanistic modelling.
  Rainfall is scaled by corresponding streamflow values raised to a power.
  This SMA uses streamflow data, so can not be used for prediction.
}
\usage{
dbm.sim(DATA, power, c = 1, qlag = 0, return_state = FALSE)
}
\arguments{
  \item{DATA}{
    time-series-like object with columns \code{P} (precipitation) and
    \code{Q} (streamflow).
  }
  \item{power}{
    power to apply to streamflow values.
  }
  \item{c}{
    constant multiplier of the result, for mass balance.
  }
  \item{qlag}{
    number of time steps to lag the streamflow (relative to rainfall)
    before multiplication.
  }
  \item{return_state}{ignored.}
}
\value{
  the simulated effective rainfall, a time series of the same length as
  the input series.  
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}(sma = "dbm")} to work with models as objects (recommended).
}
\examples{
## view default parameter ranges:
str(hydromad.options("dbm"))

data(HydroTestData)
mod0 <- hydromad(HydroTestData, sma = "dbm", routing = "expuh")
mod0

## simulate with some arbitrary parameter values
testQ <- predict(update(mod0, power = 0.5, c = 0.01, qlag = 0, tau_s = 10))
xyplot(ts.union(HydroTestData, dbm.Q = testQ))
}
\keyword{models}
