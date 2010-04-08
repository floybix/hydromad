\name{eventseq}
\alias{eventseq}
\alias{eventapply}
\alias{eventinfo}
%% old names, delete eventually:
\alias{event.clusters}
\alias{eventAttributes}
\title{
  Identify discrete events from time series and apply functions to them.
}
\description{
  Identify discrete events from time series and apply functions to them.
}
\usage{
eventseq(x, thresh = 0, inter = 1, mindur = 1, below = FALSE)

eventapply(X, events = eventseq(X, thresh, inter, mindur, below),
           FUN = sum, ..., by.column = TRUE, simplify = TRUE,
           TIMING = c("start", "middle", "end"),
           thresh = 0, inter = 1, mindur = 1, below = FALSE)

eventinfo(X, events = eventseq(X, thresh, inter, mindur, below),
          FUN = mean, ...,
          thresh = 0, inter = 1, mindur = 1, below = FALSE)
}
\arguments{
  \item{x, X}{
    a \code{\link{ts}} or \code{\link{zoo}} object.
    May be multivariate, i.e. have multiple columns.
  }
  \item{thresh}{
    
  }
  \item{inter}{
    
  }
  \item{mindur}{
    
  }
  \item{below}{

  }
  \item{events}{
    
  }
  \item{FUN, \dots}{
    
  }
  \item{by.column}{

  }
  \item{simplify}{
    if \code{FALSE}, the result will be returned as a list with one
    (named) element for each event, rather than a time series like
    object. This case allows \code{FUN} to return a complex object or
    vectors of variable lengths.
  }
  \item{TIMING}{
    
  }
}
\details{
  
}
\value{
  \code{eventseq} returns a \code{zoo} object, of the same length as
  \code{x}, with integer values indexing the identified events. Periods
  between events are \code{NA}. It has a \code{levels} attribute, so you
  can call \code{\link{nlevels}} on it to find the number of events.

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
  \code{\link[zoo]{rollapply}},
  \code{\link[latticeExtra]{panel.xblocks}},
  \code{clusters} in the \pkg{evd} package.
}
\examples{
data(SalmonBrook)
x <- window(SalmonBrook, start = "1990-01-01", end = "1995-01-01")
ev <- eventseq(x$Q, thresh = 1, inter = 5)
str(ev)
summary(ev)
eventapply(x$Q, ev, FUN = sum)
eventapply(x, ev, FUN = mean)
eventinfo(x$Q, ev)

xyplot(x) + layer_(panel.xblocks(ev, col = "grey"))
}
\keyword{ ts }
\keyword{ utilities }
