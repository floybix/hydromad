\name{IHACRES.CMD.model}
\alias{IHACRES.CMD.model}
\alias{cmd}
\alias{cmd.sim}
\title{ IHACRES Catchment Moisture Deficit (CMD) model }
\description{
  The Catchment Moisture Deficit (CMD) effective rainfall model for IHACRES.
  It is a conceptual-type model, where input rainfall is partitioned explicitly into
  drainage, evapo-transpiration, and changes in catchment moisture.
}
\usage{
cmd.sim(DATA, d, f, e, M_0 = d/2,
        return_state = FALSE)
}
\arguments{
  \item{DATA}{
    a \code{\link{ts}}-like object with named columns:
    \describe{
      \item{\code{P}}{ time series of areal rainfall depths, usually in mm. }
      \item{\code{E}}{ time series of potential evapo-transpiration, or more typically,
	temperature as an indicator of this. }
    }
  }
  \item{d}{ CMD threshold for producing flow. }
  \item{f}{ CMD stress threshold as a proportion of \code{d}. }
  \item{e}{ temperature to PET conversion factor. }
  \item{M_0}{ starting CMD value. }
  \item{return_state}{ to return state variables as well as the effective rainfall. }
}
\details{
	The mass balance step is:
	\deqn{M[t] = M[t-1] - P[t] + E_T[t] + U[t]}

	where \eqn{M} represents catchment moisture deficit (CMD),
	constrained below by 0 (the nominal fully saturated level).
	P is catchment areal rainfall, \eqn{E_T} is evapo-transpiration, and
	U is drainage (effective rainfall). All are, typically, in units of
	mm per time step.

	Rainfall effectiveness (i.e. drainage proportion) is
	a simple \emph{instantaneous} function of the CMD, with a threshold at \eqn{M = d}:
	\deqn{\frac{\mathrm{d}U}{\mathrm{d}P} = 1 - \min(1, M/d)}{
	      dU/dP = 1 - min(1, M/d)}

	The actual drainage each time step involves the integral of that relation.

	Evapo-transpiration (as a proportion of the potential rate, \code{E[t]})
	is also a simple function of the CMD, with a threshold at \eqn{M = f d}{M = f * d}:
	\deqn{E_T[t] = e E[t] \min(1, \exp\left(2\left(1 - \frac{M_f}{fd}\right)\right))}{
	      E_T[t] = e E[t] \min(1, \exp(2(1 - M_f / (fd))))}

        Note that the evapo-transpiration calculation is based on \eqn{M_f},
	which is the CMD after precipitation and drainage have been accounted for.
}
\value{
  \code{cmd.sim} returns the modelled time series of effective rainfall,
  or if \code{return_state = TRUE}, a multi-variate time series with named
  columns \code{U} (effective rainfall), \code{CMD} and
  \code{ET} (evapo-transpiration \eqn{E_T}).
}
\note{
  Normally compiled C code is used for simulation, but if
  \code{return_state = TRUE} a slower implementation in R is used.
}
\references{
Croke, B.F.W. and A.J. Jakeman (2004),
A Catchment Moisture Deficit module for the IHACRES rainfall-runoff model,
\emph{Environmental Modelling and Software}, 19(1): 1-5.

Croke, B.F.W. and A.J. Jakeman (2005),
Corrigendum to ``A Catchment Moisture Deficit module for the IHACRES rainfall-runoff model''
[Environ. Model. Softw. 19 (1) (2004) 1-5],
\emph{Environmental Modelling and Software}, 20(7): 977.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}(sma = "cmd")} to work with models as objects (recommended).
}
\examples{
## view default parameter ranges:
str(hydromad.options("cmd"))

data(Canning)
x <- cmd.sim(Canning[1:1000,], d = 200, f = 0.7, e = 0.166,
             return_state = TRUE)
xyplot(x)

data(HydroTestData)
mod0 <- hydromad(HydroTestData, sma = "cmd", routing = "expuh")
mod0

## simulate with some arbitrary parameter values
mod1 <- update(mod0, d = 200, f = 0.5, e = 0.1, tau_s = 10)
testQ <- predict(mod1, return_state = TRUE)
xyplot(cbind(HydroTestData[,1:2], cmd = testQ))
}
\keyword{ models }

