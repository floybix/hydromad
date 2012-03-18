\name{gr4j}
\encoding{UTF-8}
\alias{gr4j}
\alias{gr4j.sim}
\alias{gr4jrouting}
\alias{gr4jrouting.sim}
\title{ GR4J rainfall runoff model }
\description{
  GR4J model (mode`le du Ge´nie Rural a` 4 parame`tres Journalier).
}
\usage{
gr4j.sim(DATA, x1, etmult = 1, S_0 = 0.5, return_state = FALSE,transformed=FALSE)

gr4jrouting.sim(U, x2, x3, x4, R_0 = 0, split = 0.9,
                return_components = FALSE,
                epsilon = hydromad.getOption("sim.epsilon"),
		transformed=FALSE)
}
\arguments{
  \item{DATA}{time-series-like object with columns P (precipitation,
    mm) and E (potential evapo-transpiration, mm).}
  \item{U}{effective rainfall series.}
  \item{x1}{ maximum capacity of the production store (mm). }
  \item{x2}{ groundwater exchange coefficient (mm). }
  \item{x3}{ one day ahead maximum capacity of the routing store (mm). }
  \item{x4}{ time base of unit hydrograph UH1 (time steps). }
  \item{etmult}{Multiplier for the \code{E} input data.}
  \item{S_0}{Initial soil moisture level as fraction of \code{x1}.}
  \item{R_0}{Initial groundwater reservoir level as fraction of
\code{x3}.}
  \item{split}{ Fraction to go into quick flow routing, usually fixed at 0.9. }
  \item{return_state}{
    to return the series U, S (storage) and ET (evapotranspiration).
  }
  \item{return_components}{
    to return the series Xr, Xd and R (reservoir level).
  }
  \item{epsilon}{
    values smaller than this in the output will be set to zero.
  }
  \item{transformed}{transform parameters before use to improve identifiability. They can be untransformed using \code{\link{gr4j.transformpar}}}
}
\details{
  The default parameter ranges were taken from the "80% confidence
  intervals" given in Perrin et. al. (2003).
}
\value{
  the simulated effective rainfall, a time series of the same length as
  the input series.  
}
\references{
Perrin, C., C. Michel, et al. (2003). "Improvement of a parsimonious
model for streamflow simulation." \emph{Journal of Hydrology} 279(1-4):
275-289

\url{http://www.cemagref.fr/webgr/Modelesgb/gr4j/fonctionnement_gr4jgb.htm}
}
\author{ Felix Andrews \email{felix@nfrac.org} and Joseph Guillaume \email{josephguillaume@gmail.com} }
\seealso{
  \code{\link{hydromad}(sma = "gr4j", routing = "gr4jrouting")} to work
  with models as objects (recommended). 
}
\examples{
## view default parameter ranges:
str(c(hydromad.getOption("gr4j"),
      hydromad.getOption("gr4jrouting")))

data(HydroTestData)
mod0 <- hydromad(HydroTestData, sma = "gr4j", routing = "gr4jrouting")
mod0

## example from
## http://www.cemagref.fr/webgr/Scilab/CONT_EN/HELP_HYDROGR/c_GR4J.htm
dat <-
  cbind(P = c(0,0,0,0,0,0.04,0.59,0.03,0.01,0.16,0.37,8.76,2.65,0.05,0.02,0.02,0.38,0.00,0.02,0.46,4.46,7.71,5.71,0.79,1.33),
        E = c(0,0,0,0,0,0.24,0.24,0.24,0.24,0.24,0.25,0.25,0.26,0.27,0.28,0.32,0.33,0.34,0.35,0.36,0.36,0.37,0.37,0.38,0.38))
datz <- zoo(dat, as.Date("2000-01-01") + 1:nrow(dat))
modz <- hydromad(datz, sma = "gr4j", routing = "gr4jrouting",
    x1 = 665, x2 = 1.18, x3 = 90, x4 = 3.8, S_0 = 0.6, R_0 = 0.7)
xyplot(predict(modz, return_state = TRUE, return_components = TRUE),
       strip = FALSE, strip.left = TRUE)

## simulate with some arbitrary parameter values
mod1 <- update(mod0, x1 = 100, x2 = 20, x3 = 1, x4 = 10)
## plot results with state variables
testQ <- predict(mod1, return_state = TRUE)
xyplot(cbind(HydroTestData[,1:2], gr4j = testQ))

## show effect of increase/decrease in each parameter
parRanges <- c(hydromad.getOption("gr4j")[1],
               hydromad.getOption("gr4jrouting"))
parsims <- mapply(val = parRanges, nm = names(parRanges),
  FUN = function(val, nm) {
    lopar <- min(val)
    hipar <- max(val)
    names(lopar) <- names(hipar) <- nm
    fitted(runlist(decrease = update(mod1, newpars = lopar),
                   increase = update(mod1, newpars = hipar)))
  }, SIMPLIFY = FALSE)

xyplot.list(parsims, superpose = TRUE, layout = c(1,NA),
            strip = FALSE, strip.left = TRUE,
            main = "Simple parameter perturbation example") +
  layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
}
\keyword{models}
