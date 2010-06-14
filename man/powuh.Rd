\name{powuh}
\alias{powuh}
\alias{powuh.sim}
\alias{ssg.powuh}
\alias{normalise.powuh}
\title{ Power law transfer function models }
\description{
  Power law transfer functions
  (\emph{instantaneous unit hydrograph} or \emph{flow routing} in hydrology),
  with a single input and single output series.
}
\usage{
powuh.sim(U, delay = 0,
          a, b = 1, c = 1, init = 0, uhsteps = 100,
	  na.action = na.pass,
          epsilon = hydromad.getOption("sim.epsilon"))
}
\arguments{
  \item{U}{
    input time series.
  }
  \item{delay}{
    lag (dead time) between input and response, in time steps.
  }
  \item{a}{
    time to drop to half.
  }
  \item{b}{

  }
  \item{c}{
    curvature at half-peak point.
  }
  \item{init}{
    
  }
  \item{uhsteps}{
    number of time steps to use in approximating the unit hydrograph
    convolution filter.
  }
  \item{na.action}{
    function to remove missing values,
    e.g. \code{\link[=na.omit.ts]{na.omit}}.
  }
  \item{epsilon}{
    values smaller than this will be set to zero.
  }
}
\details{

  TODO: define model

  \code{uh <- (1 / (1 + (t/a)^(b/c))) ^ c}
  
}
\value{
  the model output as a \code{\link{ts}} object,
  with the same dimensions and time window as the input \code{U}.
}
\references{
  Croke 2006 AWR
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{expuh}}
  \code{\link{armax}}
}
\examples{

}
\keyword{ ts }

