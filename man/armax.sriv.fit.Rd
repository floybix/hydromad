\name{armax.sriv.fit}
\alias{armax.sriv.fit}
\alias{expuh.sriv.fit}
\alias{tf.sriv.fit}
\alias{tf.riv.fit}
\title{ Estimate transfer function models by Simple Refined Instrumental
  Variables method. }
\description{
  Calibrate unit hydrograph transfer function models (\code{\link{armax}}
  or \code{\link{expuh}}) using Simple Refined Instrumental
  Variables (SRIV) method.
}
\usage{
armax.sriv.fit(DATA,
           order = hydromad.getOption("order"),
           delay = hydromad.getOption("delay"),
           noise.order = hydromad.getOption("riv.noise.order"),
           fixed.ar = NULL,
           ...,
           fallback = TRUE,
           na.action = na.pass,
           epsilon = hydromad.getOption("sriv.epsilon"),
           max.iterations = hydromad.getOption("sriv.iterations"))

expuh.sriv.fit(DATA,
              order = hydromad.getOption("order"),
              delay = hydromad.getOption("delay"),
              ...,
              method = hydromad.getOption("optim.method"),
              control = hydromad.getOption("optim.control"),
              hessian = TRUE)
}
\arguments{
  \item{DATA}{
    
  }
  \item{order}{
    the transfer function order. See \code{\link{armax}}.
  }
  \item{delay}{
    delay (lag time / dead time) in number of time steps. If missing,
    this will be estimated from the cross correlation function.
  }
  \item{noise.order}{ }
  \item{fixed.ar}{ }
  \item{\dots}{
    further arguments may include
    \describe{
      \item{prefilter}{ }
      \item{initX}{ }
      \item{trace}{ ~~Describe \code{trace} here~~ }
    }
  }
  \item{fallback}{ }
  \item{na.action}{ }
  \item{epsilon}{ }
  \item{max.iterations}{ }
  \item{method, control}{
    optimisation method and settings for fitting poles constrained to
    have real, positive roots. See \code{\link{optim}}.
  }
  \item{hessian}{
    
  }
}
\details{
  In normal usage, one would not call these functions directly, but
  rather specify the routing fitting method for a \code{\link{hydromad}} model
  using that function's \code{rfit} argument. E.g. to specify fitting an
  \code{expuh} routing model by least squares one could write
  
  \code{hydromad(..., routing = "expuh", rfit = "sriv")}

  or
  
  \code{hydromad(..., routing = "expuh", rfit = list("sriv", order = c(2,1)))}.
}
\value{
  a \code{tf} object, which is a list with components
  \item{coefficients}{ the fitted parameter values.}
  \item{fitted.values}{ the fitted values. }
  \item{residuals}{ the residuals. }
  \item{delay}{ the (possibly fitted) delay time. }
}
\references{
Jakeman
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{armax}},
  \code{\link{expuh}}
}
\examples{

}
\keyword{ ts }

