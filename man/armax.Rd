\name{armax}
\alias{armax}
\alias{armax.sim}
\title{ ARMAX Transfer Function models }
\description{
	Linear transfer functions
	(\emph{instantaneous unit hydrograph} or \emph{flow routing} in hydrology),
	with a single input and single output series.
}
\usage{
armax.sim(U, a_1 = 0, a_2 = 0, a_3 = 0,
          b_0 = 1, b_1 = 0, b_2 = 0, b_3 = 0,
          pars = NULL,
	  delay = 0, init = 0, na.action = na.pass,
          epsilon = hydromad.getOption("sim.epsilon"),
          return_components = FALSE)


tfParsConvert(pars, form = c("a,b", "tau,v", "alpha,beta"))

}
\arguments{
  \item{U}{
    a \code{\link{ts}}-like object or list with named components:
    \describe{
      \item{\code{U}}{ input time series.
	For IHACRES this is effective rainfall depths, usually in mm. }
      \item{\code{Q}}{ observed reponse time series.
	For IHACRES this is discharge (streamflow) at the catchment outlet.
	It should usually be in units of mm (averaged over the catchment area).
	Use \code{\link{convertFlow}} to convert it. }
    }
  }
  \item{a_1, a_2, a_3, b_0, b_1, b_2, b_3}{
    ARMAX coefficients. Auto-regressive terms begin with \code{a} and
    moving average terms begin with \code{b}. See Details section.
  }
  \item{pars}{
    the ARMAX coefficients as a (named) vector. If this is given, it
    will over-ride the named parmameter arguments. Any number of terms
    can be given here, it is not limited to the named arguments.
  }
  \item{delay}{ lag (dead time) between input and response, in time steps. }
  \item{warmup}{ number of time steps to remove from start of record when calculating statistics. }
  \item{initX}{ whether to initialise the slower component using minimum of first 10
    time steps. }
  \item{na.action}{ function to remove missing values, e.g. \code{\link[=na.omit.ts]{na.omit}}. }
  \item{Xs_0}{ initial value of slower response component. }
  \item{Xq_0}{ initial value of quicker response component. }
  \item{return_components}{ whether to return all exponential component time series. }
}
\details{
	The transfer function used here, with input \var{u} and output \var{x} is:
	\deqn{x[t] = a_1 x[t-1] + \ldots + a_n x[t-n] +
		b_0 u[t-\delta] + \ldots + b_m u[t-m-\delta]}{
	      x[t] = a[1] x[t-1] + ... + a[n] x[t-n] +
		b[0] u[t-d] + ... + b[m] u[t-m-d]}

	and the \emph{order} is denoted \eqn{(n, m)}, with delay \eqn{\delta}{d}.

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
\seealso{ \code{\link{armax.sriv.fit}}, \code{\link{arima}} }
\examples{

}
\keyword{ ts }

