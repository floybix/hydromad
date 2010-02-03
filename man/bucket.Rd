\name{bucket}
\alias{bucket}
\alias{bucket.sim}
\title{ Single-bucket Soil Moisture Accounting models }
\description{
  Single-bucket Soil Moisture Accounting models.
}
\usage{
bucket.sim(DATA, Sb, fc, a.ei = 0, M = 1, a.ss = 0, S_0 = 0,
           return_state = FALSE)
}
\arguments{
  \item{DATA}{time-series-like object with columns P (precipitation,
    mm) and E (potential evapo-transpiration, mm).}
  \item{Sb}{Maximum soil water storage (mm).}
  \item{fc}{Field capacity (0 - 1).}
  \item{a.ei}{Interception coefficient (\eqn{\alpha_{ei}}).}
  \item{M}{Fraction of catchment area covered by deep rooted vegetation.}
  \item{a.ss}{Recession coefficients for subsurface flow from
    saturated zone (\eqn{\alpha_{ss}}).}
  \item{S_0}{Initial soil moisture level (mm).}
  \item{return_state}{to return the series U, S and ET (evapotranspiration).}
}
\details{
  From formulations given in Bai et. al. (2009), which were based on
  Farmer et. al. (2003). 
  
  The general mass balance structure is:
  \deqn{dS/dt = p - q(S) - e(S, Ep)}

  The default parameter ranges were also taken from Bai et. al. (2009).
}
\value{
  the simulated effective rainfall, a time series of the same length as
  the input series.  
}
\references{
  Farmer, D., M. Sivapalan, Farmer, D. (2003). Climate, soil and
  vegetation controls upon the variability of water balance in temperate
  and semiarid landscapes: downward approach to water balance
  analysis. \emph{Water Resources Research} 39(2), p 1035.
  
  Bai, Y., T. Wagener, P. Reed (2009). A top-down framework for
  watershed model evaluation and selection under
  uncertainty. \emph{Environmental Modelling and Software} 24(8),
  pp. 901-916.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}(sma = "bucket")} to work with models as objects (recommended).
}
\examples{
## view default parameter ranges:
str(hydromad.options("bucket"))

data(Cotter)
buc0 <- hydromad(Cotter[1:500], sma = "bucket")
buc0
simulate(buc0, 1)
}
\keyword{models}
