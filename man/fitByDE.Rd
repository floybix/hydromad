\name{fitByDE}
\alias{fitByDE}
\title{
  Fit a hydromad model using the DE (Differential Evolution) algorithm.
}
\description{
  Fit a hydromad model using the DE (Differential Evolution) algorithm.
}
\usage{
fitByDE(MODEL, 
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
    settings for the DE algorithm. See \code{\link{DEoptim.control}}.
  }
}
\details{

}
\value{
  the best model from those sampled, according to the given
  \code{objective} function.
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{DEoptim}},
  \code{\link{objFunVal}}
}
\examples{

}
\keyword{ optimization }
