\name{sacramento}
\alias{sacramento}
\alias{sacramento.sim}
\title{ Sacramento Soil Moisture Accounting model }
\description{
  Sacramento Soil Moisture Accounting model.
  Developed by the US National Weather Service.
}
\usage{
sacramento.sim(DATA, uztwm, uzfwm, uzk, pctim, adimp, zperc, rexp,
               lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree,
               etmult = 1, dt = 1, return_state = FALSE)
}
\arguments{
  \item{DATA}{time-series-like object with columns P (precipitation,
    mm) and (potential evapo-transpiration, mm).}
  \item{uztwm}{Upper zone tension water maximum capacity (mm).}
  \item{uzfwm}{Upper zone free water maximum capacity (mm).}
  \item{uzk}{Lateral drainage rate of upper zone free water expressed
    as a fraction of contents per day.}
  \item{pctim}{The fraction of the catchment which produces
    impervious runoff during low flow conditions.}
  \item{adimp}{The additional fraction of the catchment which
    exhibits impervious characteristics when the catchment's tension
    water requirements are met.}
  \item{zperc}{Maximum percolation (from upper zone free water into
    the lower zone) rate coefficient.}
  \item{rexp}{An exponent determining the rate of change of the
    percolation rate with changing lower zone water contents.}
  \item{lztwm}{Lower zone tension water maximum capacity (mm).}
  \item{lzfsm}{Lower zone supplemental free water maximum capacity (mm).}
  \item{lzfpm}{Lower zone primary free water maximum capacity (mm).}
  \item{lzsk}{Lateral drainage rate of lower zone supplemental free
    water expressed as a fraction of contents per day.}
  \item{lzpk}{Lateral drainage rate of lower zone primary free water
    expressed as a fraction of contents per day.}
  \item{pfree}{Direct percolation fraction from upper to lower zone
    free water (the percentage of percolated water which is available
    to the lower zone free water aquifers before all lower zone tension
    water deficiencies are satisfied).}
  \item{etmult}{Multiplier applied to \code{DATA$E} to estimate
    potential evapotranspiration.}
  \item{dt}{Length of each time step in days.}
  \item{return_state}{Not currently supported.}
}
\details{
  This description of the model is given by Burnash (1995):

  \dQuote{The moisture accounting system utilized in the Sacramento
    Catchment Model is a carefully structured representation of the
    catchment's soil moisture storage system. It is based on using
    simple approximations of many of those soil moisture processes
    which have been reported in the hydrologic literature. The authors
    have organised these approximations in a manner which would allow
    the determination of many catchment characteristics from carefully
    selected portions of the catchment's hydrologic record. Inasmuch as
    many of the catchment characteristics are related to the soil
    moisture capabilities of the catchment, an intelligent application
    of the model start with a good understanding of the three basic
    types of soil moisture which can potentially influence catchment
    runoff conditions. These soil moisture types are: (1) Hygroscopic
    Water, (2) Tension Water and (3) Free Water. }

  [...]

  \dQuote{Streamflow as computed by the Sacramento Catchment Model is
    the result of processing precipiatation through an algorithm
    representing the uppermost soil mantle identified as the upper zone
    and a deeper portion of the soil mantle or lower zone. The
    algorithm computes runoff in five basic forms. These are (1) direct
    runoff from permanant and temporary impervious areas, (2) surface
    runoff due to precipitation occurring at a rate faster than
    percolation and interflow can take place when both upper zone
    storages are full, (3) interflow resulting from the lateral
    drainage of a temporary free water storage, (4) supplemental base
    flow, and (5) primary base flow.} (Burnash, 1995)

  The default parameter ranges were taken from Blasone et. al. (2008).
}
\value{
  the simulated effective rainfall (\dQuote{total channel
    inflow}), a time series of the same length as the input series.
}
\references{
  Burnash, R.J.C (1995). The NWS River Forecast System -- Catchment Modeling.
  In: Vijay P. Singh (ed.), \emph{Computer models of watershed hydrology.}
  Revised edition, Highlands Ranch, Colo. : Water Resources Publications, c1995.
  \url{http://www.wrpllc.com/books/cmwh.html}.

  Blasone, R., J.A. Vrugt, H. Madsen, D. Rosbjerg, B.A. Robinson,
  G.A. Zyvoloski (2008). Generalized likelihood uncertainty estimation
  (GLUE) using adaptive Markov Chain Monte Carlo
  sampling. \emph{Advances in Water Resources} 31, pp. 630-648.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}(sma = "sacramento")} to work with models as objects (recommended).
}
\examples{
## view default parameter ranges:
str(hydromad.options("sacramento"))

data(HydroTestData)
mod0 <- hydromad(HydroTestData, sma = "sacramento")
mod0

## simulate with some arbitrary parameter values
set.seed(2)
mod1 <- simulate(update(mod0, etmult = 0.01), 1, sampletype = "random")[[1]]
testQ <- predict(mod1, return_state = TRUE)
xyplot(window(ts.union(HydroTestData[,1:2], sacramento = testQ), start = 100))
mod1
}
\keyword{ models }

