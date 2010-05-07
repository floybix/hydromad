\name{fitByDream}
\alias{fitByDream}
\title{
  Fit a hydromad model using the DE (Differential Evolution) algorithm.
}
\description{
  Fit a hydromad model using the DE (Differential Evolution) algorithm.
}
\usage{
fitByDream(MODEL, 
        objective = hydromad.getOption("objective"),
        control = hydromad.getOption("de.control"))
}
\arguments{
  \item{MODEL}{
    a model specification created by \code{\link{hydromad}}. It should
    not be fully specified, i.e one or more parameters should be defined
    by \emph{ranges} of values rather than exact values.
  }
  \item{objective}{
    objective function to minimise, given as a \code{formula} involving
    \code{Q}, \code{X} and/or \code{U}. See \code{\link{objFunVal}}.
  }
  \item{control}{
    settings for the DE algorithm. See \code{\link[DEoptim]{DEoptim.control}}.
  }
  \item{vcov}{
    if \code{vcov = TRUE}, the parameter variance-covariance matrix
    will be estimated from the last half of the sequences.
    It can be extract using \code{\link{vcov}}.
  }
}
\value{
  the best model from those sampled, according to the given
  \code{objective} function. The result from \code{\link[dream]{dream}} is
  stored as item \code{fit.result}.
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link[dream]{dream}},
  \code{\link{objFunVal}}
}
\examples{
if (require("dream")) {

  data(Cotter)
  x <- Cotter[1:1000]

  ## IHACRES CWI model with power law unit hydrograph
  modx <- hydromad(x, sma = "cwi", routing = "powuh")
  modx
  
  tmp <- fitByDream(modx, control = list(maxtime = 30))
  str(tmp$fit.result)

  summary(tmp$fit.result)

  plot(tmp$fit.result)
}
}
\keyword{ optimization }
