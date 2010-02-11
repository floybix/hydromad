\name{simulate.hydromad}
\alias{simulate.hydromad}
\title{
  Simulate hydromad models by parameter sampling.
}
\description{
  Run many simulations by sampling within parameter ranges.
}
\usage{
simulate.hydromad(object, nsim, seed, ...,
   method = c("all.combinations", "latin.hypercube", "random"),
   FUN = NULL, bind = FALSE)
x}
\arguments{
  \item{object}{
    a \code{hydromad} object (produced by the \code{\link{hydromad}()}
    function) that is not fully specified (i.e. so some parameter values
    are given as ranges).
  }
  \item{nsim}{
    number of parameter samples to run.
  }
  \item{seed}{
    optional random seed, for repeatability.
  }
  \item{\dots}{
    further arguments to \code{FUN}.
  }
  \item{method}{
    sampling method; see \code{\link{parameterSets}}.
  }
  \item{FUN}{
    optional function to apply to each simulated model. Typical examples
    would be \code{\link{objFunVal}}, \code{\link[=summary.hydromad]{summary}} or
    \code{\link[=predict.hydromad]{predict}}.
  }
  \item{bind}{
    to bind the result from \code{FUN} as one or more columns onto the
    matrix of parameter values, and return as a data frame.
  }
}
\details{
  none yet.
}
\value{
  a list of results, where each element is named by its parameter
  set. The result also has an attribute \code{"psets"} which gives the
  parameter values used in each simulation (as a data frame).

  If \code{bind = TRUE}, a data frame.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{parameterSets}},
  \code{\link{fitBySampling}}
}
\examples{
data(Canning)
mod0 <- hydromad(Canning[1:500,], sma = "cwi")
simulate(mod0, nsim = 5, method = "latin")

## plot the objective function surface over two parameters
mod1 <- update(mod0, routing = "uh", rfit = "ls")
sim1 <- simulate(mod1, 300, FUN = objFunVal, bind = TRUE)
levelplot(result ~ tw * f, sim1, cex = 2,
          panel = panel.levelplot.points)
}
\keyword{ iteration }
