\name{eventseq}
\alias{eventseq}
\alias{eventapply}
\alias{eventinfo}
\title{
  Identify discrete events from time series and apply functions to them.
}
\description{
  Identify discrete events from time series and apply functions to them.
}
\usage{
eventseq(x, thresh = 0, mingap = 1, mindur = 1,
         inthresh = thresh, inx = x,
         below = FALSE, all = FALSE, continue = FALSE,
         n = NA)

eventapply(X, events, FUN = sum, ...,
           by.column = TRUE, simplify = TRUE,
           TIMING = c("start", "middle", "end"))

eventinfo(X, events, FUN = mean, ...)
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
  \item{mingap}{
    the minimum number of time steps that can separate events. Any
    inter-event durations shorter than this will be subsumed into the
    surrounding event.
  }
  \item{mindur}{
    the minimum number of time steps in each event window. Any
    events whose duration is shorter than this will be extended.
  }
  \item{inthresh}{
    \code{inthresh} gives a second threshold value to define when
    events stop: after an event is initiated by exceeding \code{thresh},
    it continues until the second (lower) threshold \code{inthresh} is
    reached. 
  }
  \item{inx}{
    optionally, a different series may be given to determine event
    termination with \code{inthresh}. E.g. an input series may define
    event starts and a response series may define when they end.
  }
  \item{below}{
    reverses the definition of events, to identify events where the data
    falls below \code{thresh} (and/or \code{inthresh}). Setting this to
    \code{TRUE} is equivalent to negating \code{x} and \code{thresh}, etc.
  }
  \item{all}{
    to include the periods between events as additional levels (i.e. as
    events in themselves). These inter-events will be assigned negative
    levels, while the actual events will have positive levels. Otherwise
    -- in the default case -- these inter-event periods are left as
    \code{NA}. 
  }
  \item{continue}{
    set \code{continue = TRUE} to extend each event until the following
    event; otherwise there are gaps between events (the default
    behaviour). This has no effect if \code{all = TRUE}.
  }
  \item{n}{
    number of events to be identified: if this is given, then a value of
    \code{thresh} is estimated such that about this many
    events are returned. A warning is given if the actual number of
    events is not within 5 percent of \code{n}. Note that there may
    be multiple possible threshold levels giving the same number of
    events! 
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
  
  \code{eventseq} returns a zoo object, with core data consisting of an
  ordered \code{\link{factor}}, representing the identified events, and
  the same time index as \code{x}. Periods between events are left as
  \code{NA}, unless \code{all = TRUE} in which case they are treated as
  separate events.

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
  \code{\link{aggregate.zoo}},
  \code{\link{panel.xblocks}},
  \code{clusters} in the \pkg{evd} package.
}
\examples{
data(SalmonBrook)
x <- window(SalmonBrook, start = "1990-01-01", end = "1995-01-01")
ev <- eventseq(x$Q, thresh = 2, inthresh = 1)
ev2 <- eventseq(x$Q, thresh = 2, mingap = 5)

nlevels(ev) ## number of events
nlevels(ev2)
str(ev)
table(coredata(ev))
eventapply(x$Q, ev, FUN = sum)
eventapply(x, ev, FUN = mean)
eventinfo(x$Q, ev)

evplot <- xyplot(x) + 
  layer_(panel.xblocks(ev, col = "grey")) +
  layer(panel.xblocks(ev2, block.y = 0, vjust = 1, col = 1))

evplot

update(evplot, type = "s",
    xlim = as.Date(c("1990-07-01", "1990-08-31"))) +
  layer(panel.abline(h = 1:2, lty = 2), packets = 2)

## example of requesting a threshold giving about 'n' events
set.seed(0)
ee <- eventseq(rnorm(100), n = 10, mingap = 2)
nlevels(ee)

}
\keyword{ ts }
\keyword{ utilities }
