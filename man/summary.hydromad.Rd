\name{summary.hydromad}
\alias{summary.hydromad}
\alias{print.summary.hydromad}
\alias{print.summaryWithBreaks.hydromad}
\title{ Assess and summarise performance of Hydromad models }
\description{
  Assess and summarise performance of Hydromad models.
}
\usage{
\method{summary}{hydromad}(object,
	breaks = NULL, coerce = byDays,
        which = c("rel.bias", "r.squared", "r.sq.sqrt", "r.sq.log", "r.sq.monthly"),
        na.action = na.exclude, ...)
}
\arguments{
  \item{object}{ an object of class \code{hydromad}. }
  \item{breaks}{
    if specified, break up the time series and calculate statistics for each chunk.
    Can be a vector of cut points \emph{or} number of intervals to cut into \emph{or} an interval specification,
    one of "sec", "min", "hour", "day", "DSTday", "week", "month" or "year",
    optionally preceded by an integer and a space, or followed by "s".
    See \code{\link{cut.Date}}.
  }
  \item{coerce}{
    function to apply to a \code{ts} object to convert its index
    before applying \code{\link{cut}}. Typically \code{\link{byDays}},
    \code{byYears}, \code{bySecs} or \code{force} (for an integer
    index). 
  }
  \item{which}{
    which performance statistics to compute. See ...
  }
  \item{na.action}{
    an optional function to apply to the data to treat missing values before
    calculating statistics.
  }
  \item{\dots}{
    ignored.
  }
}
\details{
  The available named statistics are:
  \describe{
    \item{rel.bias}{
      
    }
    \item{...}{

    }
  }
}
\value{
  \code{summary} returns a list with named entries for each of the
  chosen stats (\code{which}). When \code{breaks} is given it retuns a
  data frame of class \code{"summaryWithBreaks"}.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}}
}
\examples{

}
\keyword{ methods }

