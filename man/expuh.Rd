\name{expuh}
\alias{expuh}
\alias{expuh.sim}
\alias{ssg.expuh}
\alias{normalise.expuh}
\title{ Exponential components transfer function models }
\description{
  Exponential components linear transfer functions
  (\emph{instantaneous unit hydrograph} or \emph{flow routing} in hydrology),
  with a single input and single output series.
}
\usage{
expuh.sim(U, delay = 0,
          tau_s = 0, tau_q = 0, tau_3 = 0,
          v_s = 1, v_q = NA, v_3 = 0, v_4 = 0,
          series = 0, loss = 0,
          Xs_0 = 0, Xq_0 = 0, X3_0 = 0,
          pars = NULL,
          return_components = FALSE,
	  na.action = na.pass,
          epsilon = hydromad.getOption("sim.epsilon"))

ssg.expuh(theta)
normalise.expuh(theta)
%tfParsConvert(pars, form = c("a,b", "tau,v", "alpha,beta"))
}
\arguments{
  \item{U}{
    input time series.
  }
  \item{delay}{
    lag (dead time) between input and response, in time steps.
  }
  \item{tau_s, tau_q, tau_3}{
    time constants for the exponential components.
  }
  \item{v_s, v_q, v_3}{
    fractional volumes for the exponential components.
  }
  \item{series}{
    defines the configuration of exponential components, as being in
    parallel and/or series (for second or third order models). See Details.
  }
  \item{loss}{
    a constant loss (or gain) term subtracted from the \emph{slow}
    (\code{s}) component.
  }
  \item{Xs_0, Xq_0, X3_0}{
    initial values of the exponential components.
  }
  \item{pars}{
    the parameters as a named vector. If this is given, it
    will over-ride the named parmameter arguments. 
  }
  \item{return_components}{
    whether to return all component time series.
  }
  \item{na.action}{
    function to remove missing values,
    e.g. \code{\link[=na.omit.ts]{na.omit}}.
  }
  \item{epsilon}{
    values smaller than this will be set to zero.
  }
  \item{theta}{
    the parameters as a named vector.
  }
}
\details{

  TODO: define model

  TODO: describe 'series'
  
}
\value{
	the model output as a \code{\link{ts}} object,
	with the same dimensions and time window as the input \code{U}.
	If \code{return_components = TRUE}, it will have multiple
	columns named
	\code{Xs}, \code{Xq} and, if relevant, \code{X3}.
}
\references{
	Jakeman, A.J., I.G. Littlewood, and P.G. Whitehead (1990),
	Computation of the instantaneous unit hydrograph and identifiable component flows
	with application to two small upland catchments,
	\emph{Journal of Hydrology}, 117: 275-300.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{expuh.sriv.fit}}
  \code{\link{armax}}
}
\examples{
data(HydroTestData)
mod1 <- hydromad(HydroTestData, routing = "expuh",
                 tau_s = 30, tau_q = 5, v_s = 0.5)
flowcomps <- predict(mod1, return_components = TRUE)
xyplot(ts.union(`Slow component` = flowcomps[,"Xs"],
                `Total flow` = flowcomps[,1] + flowcomps[,2]),
       superpose = TRUE) +
  layer(panel.refline(h = 0))
}
\keyword{ ts }

