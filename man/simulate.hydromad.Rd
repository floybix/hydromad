\name{simulate.hydromad}
\alias{simulate.hydromad}
\title{
  Simulate hydromad models by parameter sampling.
}
\description{
  Run many simulations by sampling within parameter ranges.
}
\usage{
\method{simulate}{hydromad}(object, nsim, seed, ...,
         sampletype = c("latin.hypercube", "random", "all.combinations"),
         FUN = NULL, objective = NULL, bind = !is.null(objective))
}
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
  \item{sampletype}{
    sampling method; see \code{\link{parameterSets}}.
  }
  \item{FUN}{
    optional function to apply to each simulated model. Typical examples
    would be \code{\link{objFunVal}}, \code{\link[=summary.hydromad]{summary}} or
    \code{\link[=predict.hydromad]{predict}}.
  }
  \item{objective}{
    an objective function (or statistic function); this is just an
    argument to be passed on to \code{FUN}, which in this case defaults
    to \code{\link{objFunVal}} to calculate the statistic value from
    each model. It is treated as a special argument because it is cached
    before the simulation run for efficiency.
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
simulate(mod0, nsim = 5, sampletype = "latin")

## plot the objective function surface over two parameters
mod1 <- update(mod0, routing = "armax", rfit = list("ls", order = c(2,1)))
sim1 <- simulate(mod1, 144, sampletype = "all", FUN = objFunVal,
          objective = ~ nseStat(Q, X, trans = sqrt))
levelplot(result ~ tw * f, sim1, cex = 2,
          panel = panel.levelplot.points,
          main = "R Squared (of sq.rt. data) over parameter space") +
layer_(panel.2dsmoother(...))
}
\keyword{ iteration }
