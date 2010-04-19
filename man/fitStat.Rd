\name{fitStat}
\alias{fitStat}
\alias{tsFitStat}
\alias{fitBias}
\title{ Coefficient of Efficiency (R Squared) and other fit stats }
\description{
  Coefficient of Efficiency (R Squared) and other fit stats.
}
\usage{
tsFitStat(obs, mod, ref = NULL, p = 2,
          trans = NULL, offset = identical(trans, log),
          na.action = na.pass,
          aggr = NULL, events = NULL)

fitStat(obs, mod, ref = NULL, p = 2,
        trans = NULL, offset = identical(trans, log))

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
    \code{NULL}, \code{ref} is taken as the mean of \code{obs} after
    applying any transformation (\code{trans}).
  }
  \item{p}{
    power to apply to absolute residuals (\code{abs(obs - mod)} and
    \code{abs(obs - ref)}.
  }
  \item{trans}{
    a function to apply to each data series before calculating the fit
    statistic. 
  }
  \item{offset}{
    a constant offset to apply to all data before calculating the fit
    statistic (and before applying \code{trans}). If equal to
    \code{TRUE}, the offset is taken as the 10\% quantile of the
    positive values of \code{obs}.  
  }
  \item{na.action}{
    a function to apply to the time series, which is expected to fill in
    or remove missing values (note, this is optional).
  }
  \code{aggr}{
    to aggregate the time series before calculating fit statistics
    (and before applying \code{trans} or \code{offset}). \code{aggr}
    should be a list of arguments to the \code{aggregate} function,
    which will typically dispatch to the \code{"ts"} method,
    \code{\link{aggregate.ts}}, or the \code{"zoo"} method,
    \code{\link{aggregate.zoo}}. The former accepts a \code{ndeltat}
    argument to specify the size of blocks to aggregate, and the latter
    accepts a \code{"by"} argument to specify the blocks. Common
    examples for \code{zoo} objects are \code{list(by = as.yearmon)}
    (aggregate into months), \code{list(by = cut, breaks = "3 months")}
    (aggregate into 3 months), etc.
  }
  \code{events}{
    to aggregate the time series into discrete events before calculating
    fit statistics (and before applying \code{trans} or
    \code{offset}). \code{events} should be a list of argument to the
    \code{eventseq} function. This can not be combined with \code{aggr}.
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

