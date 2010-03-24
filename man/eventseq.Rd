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
           FUN = sum, ..., by.column = TRUE,
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
  \item{TIMING}{
    
  }
}
\details{
  
}
\value{
  \code{eventseq} returns a \code{zoo} object, of the same length as
  \code{x}, with integer values indexing the identified events. Periods
  between events are \code{NA}. 

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

}
\keyword{ ts }
\keyword{ utilities }
