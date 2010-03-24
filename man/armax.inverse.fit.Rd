\name{armax.inverse.fit}
\alias{armax.inverse.sim}
\alias{armax.inverse.fit}
\alias{expuh.inverse.fit}
\alias{tf.inverse.fit}
\title{ Estimate transfer function models by Inverse Filtering. }
\description{
  Calibrate unit hydrograph transfer function models (\code{\link{armax}}
  or \code{\link{expuh}}) using Inverse Filtering.
}
\usage{
armax.inverse.sim(Q, P = NULL,
                  pars = c(a_1 = 0, b_0 = 1),
                  delay = 0, Xs_0 = 0, Xq_0 = 0,
                  rain.factor = 1.1,
                  rises.only = FALSE,
                  use.Qm = TRUE,
                  use.fft.method = FALSE,
                  constrain.fft = TRUE,
                  mass.balance = use.fft.method,
                  scale.window = NA)

armax.inverse.fit(DATA,
                  order = hydromad.getOption("order"),
                  delay = hydromad.getOption("delay"), 
                  fit.method = hydromad.getOption("inverse.fit.method"),
                  normalise = TRUE, 
                  init.U = TRUE, pars = NULL,
                  use.Qm = TRUE,
                  fft.inverse.sim = FALSE, 
                  rises.only = FALSE,
                  ...,
                  max.iterations = hydromad.getOption("inverse.iterations"), 
                  rel.tolerance = hydromad.getOption("inverse.rel.tolerance"), 
                  par.epsilon = hydromad.getOption("inverse.par.epsilon"), 
                  init.attempt = 0, trace = hydromad.getOption("trace"))

expuh.inverse.fit(DATA,
                  order = hydromad.getOption("order"),
                  delay = hydromad.getOption("delay"),
                  ...)
}
\arguments{
  \item{DATA}{
    a \code{\link{ts}}-like object with named columns:
    \describe{
      \item{\code{Q}}{ observed output time series. }
      \item{\code{P}}{ observed input time series (optional). }
    }
  }
  \item{order}{
    the transfer function order. See \code{\link{armax}}.
  }
  \item{delay}{
    delay (lag time / dead time) in number of time steps. If missing,
    this will be estimated from the cross correlation function.
  }



  
}
\details{
  In normal usage, one would not call these functions directly, but
  rather specify the routing fitting method for a \code{\link{hydromad}} model
  using that function's \code{rfit} argument. E.g. to specify fitting an
  \code{expuh} routing model by inverse filtering one could write
  
  \code{hydromad(..., routing = "expuh", rfit = "inverse")}

  or
  
  \code{hydromad(..., routing = "expuh", rfit = list("inverse", order = c(2,1)))}.
}
\value{
  a \code{tf} object, which is a list with components
  \item{coefficients}{ the fitted parameter values.}
  \item{fitted.values}{ the fitted values. }
  \item{residuals}{ the residuals. }
  \item{delay}{ the (possibly fitted) delay time. }
}
\references{
...
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{armax}},
  \code{\link{expuh}},
  \code{\link{armax.sriv.fit}},
}
\examples{
U <- ts(c(0, 0, 0, 1, rep(0, 30), 1, rep(0, 20)))
Y <- expuh.sim(lag(U, -1), tau_s = 10, tau_q = 2, v_s = 0.5, v_3 = 0.1)
set.seed(0)
Yh <- pmax(Y + rnorm(Y, sd = 0.01), 0)
fit1 <- armax.inverse.fit(ts.union(P = U, Q = Yh),
                          order = c(2, 2), warmup = 0)
fit1
xyplot(ts.union(observed = Yh, fitted = fitted(fit1)),
       superpose = TRUE)
}
\keyword{ ts }

