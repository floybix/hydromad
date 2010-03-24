\name{lambda}
\alias{lambda}
\alias{lambda.sim}
\alias{ssg.lambda}
\alias{normalise.lambda}
\title{ Two exponential components transfer function with variable partitioning }
\description{
  Two exponential components transfer function with variable partitioning.
}
\usage{
lambda.sim(U, delay = 0,
           tau_s = 0, tau_q = 0,
           lambda = 0, v_s = 1,
           loss = 0,
           Xs_0 = 0, Xq_0 = 0,
           return_components = FALSE,
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
  \item{tau_s, tau_q}{
    time constants for the exponential components.
  }
  \item{lambda}{
    
  }
  \item{v_s}{
    maximum fractional volume for the slower exponential component.
  }
  \item{loss}{
    a constant loss (or gain) term subtracted from the \emph{slow}
    (\code{s}) component.
  }
  \item{Xs_0, Xq_0}{
    initial values of the exponential components.
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
}
\details{

  TODO: define model

  
}
\value{
  the model output as a \code{\link{ts}} object,
  with the same dimensions and time window as the input \code{U}.
  If \code{return_components = TRUE}, it will have multiple
  columns named \code{Xs} and \code{Xq}.
}
\references{
  ...
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{expuh}}
  \code{\link{lambda.inverse.fit}}
}
\examples{

}
\keyword{ ts }

