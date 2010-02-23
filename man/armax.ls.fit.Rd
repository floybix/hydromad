\name{armax.ls.fit}
\alias{armax.ls.fit}
\alias{expuh.ls.fit}
\alias{tf.ls.fit}
\title{ Estimate transfer function models by Least Squares. }
\description{
  Calibrate unit hydrograph transfer function models (\code{\link{armax}}
  or \code{\link{expuh}}) using Least Squares with prefiltering.
}
\usage{
armax.ls.fit(DATA,
             order = hydromad.getOption("order"),
             delay = hydromad.getOption("delay"),
             prefilter = hydromad.getOption("prefilter"),
             warmup = hydromad.getOption("warmup"),
             normalise = FALSE,
             fixed.ar = NULL,
             weights = NULL,
             initX = TRUE,
             na.action = na.pass,
             trace = hydromad.getOption("trace"))

expuh.ls.fit(DATA,
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
  \item{prefilter}{
  }
  \item{warmup}{
  }
  \item{normalise}{
  }
  \item{fixed.ar}{
  }
  \item{weights}{
  }
  \item{initX}{
  }
  \item{na.action}{
  }
  \item{trace}{
  }
  \item{\dots}{
    passed on to \code{armax.ls.fit}.
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
  
  \code{hydromad(..., routing = "expuh", rfit = "ls")}

  or
  
  \code{hydromad(..., routing = "expuh", rfit = list("ls", order = c(2,1)))}.
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
  \code{\link{expuh}},
  \code{\link{armax.sriv.fit}},
  \code{\link{arima}}
}
\examples{

}
\keyword{ ts }

