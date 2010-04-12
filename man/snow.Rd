\name{snow}
\alias{snow}
\alias{snow.sim}
\title{ Simple degree day factor snow model }
\description{
  Simple degree day factor snow model coupled with IHACRES CMD soil
  moisture model. 
}
\usage{
snow.sim(DATA, Tmax, Tmin, kd, kf, rcap, Tmelt = Tmin,
         cr = 1, cs = 1, LSWE_0 = 0, ISWE_0 = 0,
         d, f, e, M_0 = d/2, return_state = FALSE)
}
\arguments{
  \item{DATA}{ a \code{\link{ts}}-like object with named columns:
    \describe{
      \item{\code{P}}{ time series of areal rainfall depths, usually in mm. }
      \item{\code{E}}{ time series of potential evapo-transpiration, or more typically,
	temperature as an indicator of this. }
    }
  }
  \item{Tmax}{
    temperature threshold for rain, all rain is liquid above this
    threshold.
  }
  \item{Tmin}{
    temperature threshold for rain, all rain is snow below this
    threshold.
  }
  \item{Tmelt}{
    temperature threshold for snowmelt and freezing in the snowpack.
  }
  \item{kd}{
    degree day factor for snowmelt.
  }
  \item{kf}{
    degree day factor for freezing.
  }
  \item{rcap}{
    retention parameter for liquid water capacity of snowpack.
  }
  \item{cr}{ correction factor for rainfall. }
  \item{cs}{ correction factor for snowfall. }
  \item{d, f, e, M_0}{
    parameters for the \link{IHACRES.CMD.model}.
  }
  \item{return_state}{
    to return state variables as well as the effective rainfall.
  }
}
\details{
  SWE snow water equivalent
  
  ISWE water equivalent of ice in the snowpack
  
  LSWE liquid water retained in the snowpack

}
\value{
  \code{snow.sim} returns the modelled time series of effective rainfall,
  or if \code{return_state = TRUE}, a multi-variate time series with named
  columns \code{U} (effective rainfall), \code{SWE} (snow water
  equivalent) and \code{TF}, as well as the CMD state variables.
}
\references{
Kokkonen T., Jakeman A.J, Koivusalo.H, Norton.J.:
 COMPUTATIONAL METHODS FOR WATER RESOURCE ASSESSMENTS:
 AN EXERCISE KIT
 Educational Series on Modelling and Software
 iEMSs International Modelling and Software Society
 Available through www.iemss.org
}
\author{
  Coded in R by Jarkko Koskela @tkk.fi 2010-02-26.

  Converted to C by Felix Andrews \email{felix@nfrac.org}.
}
\seealso{
  \code{\link{hydromad}(sma = "snow")} to work with models as objects (recommended).
}
\examples{
## view default parameter ranges:
str(hydromad.options("snow"))

data(HydroTestData)
mod0 <- hydromad(HydroTestData, sma = "snow", routing = "expuh")
mod0

## simulate with some arbitrary parameter values
mod1 <- update(mod0, Tmax = 15, Tmin = 5, cr = 1, cs = 1, 
               kd = 2, kf = 0, rcap = 0.5,
               d = 200, f = 0.5, e = 0.1, tau_s = 10)
testQ <- predict(mod1, return_state = TRUE)
xyplot(cbind(HydroTestData[,1:2], snow = testQ))
}
\keyword{ models }

