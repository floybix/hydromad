\name{fitStat}
\alias{fitStat}
\alias{fitBias}
\title{ Coefficient of Efficiency (R Squared) and other fit stats }
\description{
  Coefficient of Efficiency (R Squared) and other fit stats.
}
\usage{
fitStat(obs, mod, ref = NULL, trans = NULL,
        offset = identical(trans, log), p = 2)

fitBias(obs, mod, relative.bias = TRUE, na.rm = TRUE)
}
\arguments{
  \item{obs}{
    observed data vector.
  }
  \item{mod}{
    model-predicted data vector corresponding to \code{obs}.
  }
  \item{ref}{
    reference model predictions corresponding to \code{obs}. If
    \code{NULL}, \code{ref} is taken as a constant value: the mean of
    \code{obs}. 
  }
  \item{trans}{
    a function to apply to all data before calculating the fit
    statistic. 
  }
  \item{offset}{
    a constant offset to apply to all data before calculating the fit
    statistic. If equal to \code{TRUE}, the offset is taken as the 10\%
    quantile of the positive values of \code{obs}. 
  }
  \item{p}{
    power to apply to absolute residuals (\code{abs(obs - mod)} and
    \code{abs(obs - ref)}.
  }
  \item{relative.bias}{
    if \code{FALSE}, the result is the mean error in the original data
    units: \code{mean(mod - obs)}. Otherwise, the result is divided by
    \code{mean(obs)} to give the relative bias.
  }
  \item{na.rm}{
    to ignore missing values.
  }
}
\details{
  
}
\value{
  a single numeric value.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{summary.hydromad}},
  \code{\link{objFunVal}}
}
\examples{

}
\keyword{ ts }

