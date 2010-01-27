\name{simulate.ihacres}
\alias{simulate.ihacres}
\title{
  Simulate ihacres models by parameter sampling.
}
\description{
  Run many simulations by sampling within parameter ranges.
}
\usage{
simulate.ihacres(object, nsim, seed, ...,
   method = c("all.combinations", "latin.hypercube", "random"),
   FUN = NULL)
}
\arguments{
  \item{object}{
    an \code{ihacres} object (produced by the \code{\link{ihacres}()}
    function) that is not fully specified (i.e. so some parameter values
    are given as ranges).
  }
  \item{nsim}{
    number of parameter samples to run.
  }
  \item{seed}{
    optional random seed, for repeatability.
  }
  \item{method}{
    sampling method; see \code{\link{parameterSets}}.
  }
  \item{FUN}{
    optional function to apply to each simulated model. Typical examples
    would be \code{\link[summary.ihacres]{summary}} or
    \code{\link[predict.ihacres]{predict}}.
}
  \item{\dots}{
    further arguments to \code{FUN}.
  }
}
\details{
  none yet.
}
\value{
  a list of results, where each element is named by its parameter
  set. The result also has an attribute \code{"psets"} which gives the
  parameter values used in each simulation (as a data frame).
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{parameterSets}},
  \code{\link{fitBySampling}}
}
\examples{
data(Canning)
x <- Canning["1980"]
mod0 <- ihacres(x, sma = "cwi")
simulate(mod0, nsim = 12, method = "latin")
}
\keyword{ iteration }
