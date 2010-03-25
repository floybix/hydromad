\name{runoffratio}
\alias{runoffratio}
\alias{runoffratio.sim}
\alias{absorbScale.hydromad.runoffratio}
\title{ Simple time-varying runoff proportion }
\description{
  Simple time-varying runoff proportion. Rainfall is scaled by the
  runoff coefficient estimated in a moving window.
  This SMA uses streamflow data, so can not be used for prediction.
}
\usage{
runoffratio.sim(DATA, width = 30, kernel = 2, sides = 2, rrthresh = 0,
                qlag = 0, scale = 1, return_state = FALSE)
}
\arguments{
  \item{DATA}{
    time-series-like object with columns \code{P} (precipitation) and
    \code{Q} (streamflow).
  }
  \item{width}{
    width of the time window (in time steps) in which to estimate the
    runoff coefficient. 
  }
  \item{kernel}{
    type of window used to estimate the runoff coefficient: 1 is
    rectangular, 2 is triangular-weighted, 3 is Gaussian-like. 
  }
  \item{sides}{
    2 for time-centered estimates, 1 for estimates using data backward
    in time only.
  }
  \item{rrthesh}{
    a theshold value of the runoff ratio, below which there is no
    effective rainfall.
  }
  \item{qlag}{
    number of time steps to lag the streamflow (relative to rainfall)
    before estimating the runoff coefficient.
  }
  \item{scale}{
    constant multiplier of the result, for mass balance.
    If this parameter is set to \code{NA} (as it is by default)
    in \code{\link{hydromad}} it will be set by mass balance calculation. 
  }
  \item{return_state}{ignored.}
}
\value{
  the simulated effective rainfall, a time series of the same length as
  the input series.  
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}(sma = "runoffratio")} to work with models as objects (recommended).
}
\examples{
## view default parameter ranges:
str(hydromad.options("runoffratio"))

data(HydroTestData)
mod0 <- hydromad(HydroTestData, sma = "runoffratio", routing = "expuh")
mod0

## simulate with some arbitrary parameter values
testQ <- predict(update(mod0, width = 30, rrthresh = 0.2, tau_s = 10))
xyplot(ts.union(HydroTestData, runoffratio.Q = testQ))
}
\keyword{models}
