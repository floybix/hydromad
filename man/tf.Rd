\name{tf}
\alias{tf}
\alias{xferfun}
\alias{tf.sim}
\alias{predict.tf}
\alias{update.tf}
\alias{coef.tf}
\alias{summary.tf}
\alias{print.summary.tf}
\alias{print.tf}
\alias{normalise.tf}
\alias{vcov.tf}
\title{ Transfer Function models }
\description{
	Linear transfer functions
	(\emph{instantaneous unit hydrograph} or \emph{flow routing} in hydrology),
	with a single input and single output series.
}
\usage{
tf(DATA = list(U=, Q=),
	pars = c(tau_s = , tau_q = 0, v_s = 1, v_3 = 0),
	delay = 0, warmup, initX = TRUE,
	na.action = na.pass)

tf.sim(U,
	pars = c(tau_s = , tau_q = 0, v_s = 1, v_3 = 0),
	delay = 0, Xs_0 = 0, Xq_0 = 0,
	return_components = FALSE,
	na.action = na.pass, epsilon)

% move into tf.object.Rd ?
% TODO: note methods for coef, residuals, fitted, vcov
\method{summary}{tf}(object,
	which = c("arpe", "rel.bias", "abs.bias",
        "r.squared", "r.sq.sqrt", "r.sq.log", "r.sq.monthly",
	"ssg", "residuals"), ...)

\method{coef}{tf}(object, form = c("a,b", "tau,v", "alpha,beta"), ...)

\method{predict}{tf}(object, newdata = NULL,
	return_components = FALSE, ...)

\method{update}{tf}(object, newdata = NULL, pars, ...)

normalise.tf(object, ...)
}
\arguments{
  \item{DATA, U}{ a \code{\link{ts}}-like object or list with named components:
	\describe{
	  \item{\code{U}}{ input time series.
	  	For IHACRES this is effective rainfall depths, usually in mm. }
	  \item{\code{Q}}{ observed reponse time series.
		For IHACRES this is discharge (streamflow) at the catchment outlet.
		It should usually be in units of mm (averaged over the catchment area).
		Use \code{\link{convertFlow}} to convert it. }
	    }
	    \emph{DATA} may be omitted, to create an empty model object. }
  \item{pars}{ specification of the model; see Details section. }
  \item{delay}{ lag (dead time) between input and response, in time steps. }
  \item{warmup}{ number of time steps to remove from start of record when calculating statistics. }
  \item{initX}{ whether to initialise the slower component using minimum of first 10
    time steps. }
  \item{na.action}{ function to remove missing values, e.g. \code{\link[=na.omit.ts]{na.omit}}. }
  \item{Xs_0}{ initial value of slower response component. }
  \item{Xq_0}{ initial value of quicker response component. }
  \item{return_components}{ whether to return all exponential component time series. }
  \item{object}{ }
  \item{which}{ }
  \item{newdata}{ }
  \item{\dots}{ }
}
\details{
	The transfer function used here, with input \var{u} and output \var{x} is:
	\deqn{x[t] = a_1 x[t-1] + \ldots + a_n x[t-n] +
		b_0 u[t-\delta] + \ldots + b_m u[t-m-\delta]}{
	      x[t] = a[1] x[t-1] + ... + a[n] x[t-n] +
		b[0] u[t-d] + ... + b[m] u[t-m-d]}

	and the \emph{order} is denoted \eqn{(n, m)}, with delay \eqn{\delta}{d}.

%	\deqn{x[t] = \sum_{i=1}^n \frac{\beta_i}{1 + \alpha_i z^{-1}} u[t]}
%	\eqn{z^{-1}} is the backshift operator.

	In most cases this can be expressed as a sum of exponentially decaying components,
	each defined by a recession rate \eqn{\alpha} and peak response \eqn{\beta}.
	If there are two components, which we call quick \eqn{x_q} and slow \eqn{x_s}, then
	\eqn{x[t] = x_q[t] + x_s[t]}, and:
	\deqn{x_q[t] = \alpha_q x_q[t-1] + \beta_q u[t]}{x_q[t] = \alpha_q . x_q[t-1] + \beta_q . u[t]}
	\deqn{x_s[t] = \alpha_s x_s[t-1] + \beta_s u[t]}{x_s[t] = \alpha_s . x_s[t-1] + \beta_s . u[t]}

	It is often easier to define these exponential components in terms of
	time constants \eqn{\tau} (number of time steps to reduce to a fraction
	\eqn{1/e \approx 37\%}) and fractional volumes \eqn{v}:
	\deqn{\tau = -1 / \log(\alpha)}
	\deqn{v = \beta / (1 - \alpha)}
	In this form it is assumed that the total volume (\emph{steady state gain}) is 1.

	\code{tf} and \code{tf.sim} allow the model to be specified in
	either \eqn{\{\tau, v\}} form, or \eqn{\{a, b\}} form. In the former, four model structures can
	be specified via the \code{pars} argument:
	\itemize{
		\item \code{(tau_s=)}, a single exponential store.
		\item \code{(tau_s=, v_s=)}, exponential store and instantaneous component in parallel.
			The instantaneous store has volume \code{v_q = 1 - v_s}.
		\item \code{(tau_s=, tau_q=, v_s=)}, two exponential stores in parallel.
		\item \code{(tau_s=, tau_q=, v_s=, v_3=)}, two exponential stores and an instantaneous component in parallel.
			Here \code{v_q = 1 - v_s - v_3} and \code{tau_3 = 0}.
	}

	Alternatively, models can be specified in \eqn{\{a, b\}} form, like
	\code{(a_1=, a_2=, a_3=, b_0=, b_1=, b_2=)}.
}
\value{
	\code{tf} returns an object of class \code{tf}.

	\code{tf.sim} returns the model output as a \code{\link{ts}} object,
	with the same dimensions and time window as the input \code{U}.
	If \code{return_components = TRUE}, a list is returned with components
	\code{Xs}, \code{Xq} and, if relevant, \code{X3}.
	These are all \code{ts} objects as above.
}
\references{
	Jakeman, A.J., I.G. Littlewood, and P.G. Whitehead (1990),
	Computation of the instantaneous unit hydrograph and identifiable component flows
	with application to two small upland catchments,
	\emph{Journal of Hydrology}, 117: 275-300.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{ \code{\link{tf.fit}}, \code{\link{arima}} }
\examples{

}
\keyword{ ts }

