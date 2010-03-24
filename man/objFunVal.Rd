\name{objFunVal}
\alias{objFunVal}
\title{ Calculate objective function value for a fitted model. }
\description{
  Calculate objective function value for a fitted model.
}
\usage{
objFunVal(model, objective = hydromad.getOption("objective"),
          nan.ok = FALSE)
}
\arguments{
  \item{model}{
    a fitted \code{hydromad} model object, i.e. it must have specific
    parameter values, not ranges.
  }
  \item{objective}{
    the objective function, as a \code{formula}. See Details.
  }
  \item{nan.ok}{
    by default, an error is thrown if the result is \code{NaN}. Set this
    argument to \code{TRUE} to suppress the error.
  }
}
\details{
 
  Q

  X

  
  
  
}
\value{
  the objective function value. It must be numeric and of length one;
  anything else is an error.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \link{hydromad.object},
  \code{\link{hydromad.options}}
}
\examples{

}
\keyword{ ts }

