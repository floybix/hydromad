\name{fitStat}
\alias{fitStat}
\alias{tsFitStat}
\alias{fitBias}
\title{ Coefficient of Efficiency (R Squared) and other fit stats }
\description{
  Coefficient of Efficiency (R Squared) and other fit stats.
}
\usage{
tsFitStat(obs, mod, ref = NULL, ...,
          na.action = na.pass,
          aggr = NULL, events = NULL)

fitStat(obs, mod, ref = NULL, p = 2,
        trans = NULL, offset = identical(trans, log),
        negatives.ok = FALSE)

fitBias(obs, mod, relative.bias = TRUE)
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
  \item{...}{
    passed to \code{fitStat}.
  }
  \item{na.action}{
    a function to apply to the time series, which is expected to fill in
    or remove missing values (note, this is optional).
  }
  \item{aggr}{
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
  \item{events}{
    to aggregate the time series into discrete events before calculating
    fit statistics (and before applying \code{trans} or
    \code{offset}). \code{events} should be a list of argument to the
    \code{eventseq} function. This can not be combined with \code{aggr}.
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
  \item{negatives.ok}{
    if \code{FALSE}, the default case, all values in \code{obs},
    \code{mod} and \code{ref} are constrained to be non-negative;
    i.e. negative values are replaced with zero.
  }
  \item{relative.bias}{
    to calculate the relative bias, i.e. \code{mean(mod - obs) /
      mean(obs)} (except that missing values are omitted).
  }
}
\details{
  The \code{tsFitStat} function will merge its arguments together, such
  that corresponding time steps are compared to each other even if the
  time windows are not equal. In \code{fitStat}, it is assumed that the
  series are already synchronised.
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

