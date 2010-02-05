\name{uh.riv.fit}
\alias{uh.riv.fit}
\alias{expuh.riv.fit}
\alias{tf.riv.fit}
\title{ Estimate transfer function models by Refined Instrumental
  Variables method. }
\description{
  Calibrate unit hydrograph transfer function models (\code{\link{uh}}
  or \code{\link{expuh}}) using Refined Instrumental
  Variables (RIV) or Simple Refined Instrumental
  Variables (SRIV) method.
}
\usage{
uh.riv.fit(DATA,
          order = hydromad.getOption("order"),
          delay = hydromad.getOption("delay"),
          ...)

expuh.riv.fit(DATA,
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
    the transfer function order. See \code{\link{uh}}.
  }
  \item{delay}{
    delay (lag time / dead time) in number of time steps. If missing,
    this will be estimated from the cross correlation function.
  }
  \item{\dots}{
    further arguments may include
    \describe{
      \item{prefilter}{ ~~Describe \code{prefilter} here~~ }
      \item{na.action}{ ~~Describe \code{na.action} here~~ }
      \item{trace}{ ~~Describe \code{trace} here~~ }
    }
  }
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
  
  \code{hydromad(..., routing = "expuh", rfit = "riv")}

  or
  
  \code{hydromad(..., routing = "expuh", rfit = list("riv", order = c(2,1)))}.
}
\value{
  a list with components
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
  \code{\link{uh}},
  \code{\link{expuh}},
  \code{\link{uh.ls.fit}},
  \code{\link{arima}}
}
\examples{

}
\keyword{ ts }

