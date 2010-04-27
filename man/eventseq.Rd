\name{eventseq}
\alias{eventseq}
\alias{eventapply}
\alias{eventinfo}
\alias{levels.eventseq}
\title{
  Identify discrete events from time series and apply functions to them.
}
\description{
  Identify discrete events from time series and apply functions to them.
}
\usage{
eventseq(x, thresh = 0, mingap = 1, mindur = 1,
         below = FALSE, all = FALSE, inter = NA)

eventapply(X, events, FUN = sum, ...,
           by.column = TRUE, simplify = TRUE,
           TIMING = c("start", "middle", "end"))

eventinfo(X, events, FUN = mean, ...)

\method{levels}{eventseq}(x)
}
\arguments{
  \item{x, X}{
    a \code{\link{ts}} or \code{\link{zoo}} object.
    May be multivariate, i.e. have multiple columns.
    In this case, events are defined jointly over the multiple
    time series: events continue while any series exceeds the
    threshold (or remains below the threshold if \code{below = TRUE}. 
  }
  \item{thresh}{
    threshold value: the data must be strictly above this level (or, if
    \code{below = TRUE}, below this level) to define an event. 
  }
  \item{mingap, inter}{
    the minimum number of time steps that can separate events. Any
    inter-event durations shorter than this will be subsumed into the
    surrounding event. (\code{inter} is the old name, now deprecated).
  }
  \item{mindur}{
    the minimum number of time steps required to define an event. Any
    events whose duration is shorter than this will be ignored (not
    identified as an event).
  }
  \item{below}{
    to identify events where the data falls below \code{thresh}.
  }
  \item{all}{
    to include the periods between events as additional levels (i.e. as
    events in themselves). These inter-events will be assigned negative
    levels, while the actual events will have positive levels. Otherwise
    -- in the default case -- these inter-event periods are left as
    \code{NA}. 
  }
  \item{events}{
    a \code{factor}-like object, perhaps \code{\link{ts}} or
    \code{\link{zoo}}. Its levels (unique values) define events, and
    \code{NA} values are ignored. If \code{events} is a factor or
    vector, its values are assumed to correspond to \code{X}; otherwise,
    it is assumed to be a time series object and is merged
    (\code{cbind}ed) with \code{X}.
  }
  \item{FUN, \dots, by.column}{
    a function to apply to the values in each event. In the default case
    of \code{by.column = TRUE}, \code{FUN} will be passed a vector of
    values from one column of \code{X}. When \code{by.column = FALSE},
    which only makes sense when \code{X} is a matrix, \code{FUN} will be
    passed a matrix. The dots (\code{\dots}) are passed on too.
  }
  \item{simplify}{
    if \code{FALSE}, the result will be returned as a list with one
    (named) element for each event, rather than a time series like
    object. This case allows \code{FUN} to return a complex object or
    vectors of variable lengths.
  }
  \item{TIMING}{
    defines how to construct the time index of the result. Should the
    time corresponding to an event be taken from the \code{time()} of
    its start, middle, or end?
  }
}
\value{
  
  \code{eventseq} returns a vector of the same class as its input
  \code{x}, plus \code{"eventseq"}, with integer values indexing the
  identified events. Periods between events are left as \code{NA},
  unless \code{all = TRUE} in which case they will have unique negative
  integer values. The result has a \code{levels} method, so you can
  call \code{\link{nlevels}} on it to find the number of events.

  \code{eventapply} returns a \code{zoo} object (an irregular time
  series in this case), with the value returned from \code{FUN}
  applied to each discrete event in \code{X}.

  \code{eventinfo} returns a \code{data.frame} with columns
  \describe{
    \item{\code{Time}}{
      time that the event started (from \code{time(X)}).
    }
    \item{\code{Month, Year}}{
      month and year (as integers) of the mid-point of the event.
    }
    \item{\code{Value}}{
      result of \code{FUN} applied to the event.
    }
    \item{\code{Duration}}{
      length of the event in time steps / data points.
    }
    \item{\code{PreDuration}}{
      number of time steps since the last event ended.
    }
  }
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{cut.Date}},
  \code{\link{tapply}},
  \code{\link{rollapply}},
  \code{\link{panel.xblocks}},
  \code{clusters} in the \pkg{evd} package.
}
\examples{
data(SalmonBrook)
x <- window(SalmonBrook, start = "1990-01-01", end = "1995-01-01")
ev <- eventseq(x$Q, thresh = 1, mingap = 5)
str(ev)
summary(ev)
eventapply(x$Q, ev, FUN = sum)
eventapply(x, ev, FUN = mean)
eventinfo(x$Q, ev)

xyplot(x) + layer_(panel.xblocks(ev, col = "grey"))
}
\keyword{ ts }
\keyword{ utilities }
