\name{fitByCMAES}
\alias{fitByCMAES}
\title{
  Fit a hydromad model using CMA-ES (Covariance matrix adapting
  evolutionary strategy) from cmaes package
}
\usage{
fitByCMAES(MODEL, objective = hydromad.getOption("objective"), control = hydromad.getOption("cmaes.control"), vcov = FALSE)
}
\arguments{
  \item{MODEL}{
        a model specification created by \code{\link{hydromad}}. It should
    not be fully specified, i.e one or more parameters should be defined
    by \emph{ranges} of values rather than exact values.
}
\item{objective}{
  objective function to maximise, given as a \code{function(Q, X, ...)}.
  See \code{\link{objFunVal}}.
}
  \item{control}{
    settings for the CMA-ES algorithm. See \code{\link{cma_es}}
}
  \item{vcov}{
    Ignored
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    the best model from those sampled, according to the given
  \code{objective} function. Also, these extra elements are inserted:
  \item{fit.result}{
    the result from \code{\link{SCEoptim}}.
  }
  \item{objective}{
    the \code{objective} function used.
  }
  \item{funevals}{
    total number of evaluations of the model simulation function.
  }
  \item{timing}{
    timing vector as returned by \code{system.time}.
  }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Joseph Guillaume
}
\seealso{
  \code{\link{cma_es}}
}
\keyword{optimisation}