\name{summary.hydromad}
\alias{summary.hydromad}
\alias{print.summary.hydromad}
\alias{print.summaryWithBreaks.hydromad}
\alias{summary.hydromad.runlist}
\title{ Assess and summarise performance of Hydromad models }
\description{
  Assess and summarise performance of Hydromad models.
}
\usage{
\method{summary}{hydromad}(object,
	breaks = NULL, coerce = byDays,
        stats = c("rel.bias", "r.squared", "r.sq.sqrt", "r.sq.log", "r.sq.monthly"),
        na.action = na.exclude, ...)

\method{summary}{hydromad.runlist}(object, ...,
        stats = c("rel.bias", "r.squared", "r.sq.sqrt", "r.sq.log", "r.sq.monthly"),
        items = stats)
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
    \code{byYears}, \code{bySecs} or \code{identity} (for an integer
    index). 
  }
  \item{stats}{
    which performance statistics to compute. See Details.
  }
  \item{items}{
    see \code{\link{summary.runlist}}.
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
    \item{bias, rel.bias}{
      bias in data units, \eqn{\sum ( X - Q )},
      \code{fitBias(Q, X, rel = FALSE)}.
    }
    \item{rel.bias}{
      bias in data units, \eqn{\sum ( X - Q ) / \sum Q},
      \code{fitBias(Q, X)}.     
    }
    \item{r.squared}{
      R Squared (coefficient of determination),
      \eqn{1 - \sum (Q - X)^2 / \sum (Q - \bar{Q})^2},
      \code{fitStat(Q, X)}.
    }
    \item{r.sq.sqrt}{
      R Squared (coefficient of determination) using square-root
      transformed data (less weight on peak flows),
      \eqn{1 - \frac{\sum |\sqrt{Q} - \sqrt{X}|^2 }{ \sum |\sqrt{Q} - \bar{\sqrt{Q}}|^2 }},
      \code{fitStat(Q, X, trans = sqrt)}.
    }
    \item{r.sq.log}{
      R Squared (coefficient of determination) using log
      transformed data (roughly equal weight to all flow levels),
      \eqn{1 - \frac{\sum |\log{(Q+\epsilon)} - \log{(X+\epsilon)}|^2 }{ \sum |\log{(Q+\epsilon)} - \bar{\log{(Q+\epsilon)}}|^2 }},
      \code{fitStat(Q, X, trans = log)}.
      By default \eqn{\epsilon} is the 90 percentile (i.e. lowest
      decile) of the non-zero values of Q.
    }
    \item{r.sq.monthly}{
      R Squared (coefficient of determination) with data aggregated into
      blocks of 30 time steps, 
      \code{fitStat(Q, X, aggr=30)}.
    }
    \item{U1}{

    }
    \item{X1}{
      correlation of model with error lagged one time step
      \code{cor(head(X,-1), tail(X-Q,-1))}
    }
    \item{arpe, yic}{
      
    }
  }
}
\value{
  \code{summary} returns a list with named entries for each of the
  chosen stats (\code{stats}). When \code{breaks} is given it retuns a
  \code{zoo} object of class \code{"summaryWithBreaks"}.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}}
}
\examples{

}
\keyword{ methods }

