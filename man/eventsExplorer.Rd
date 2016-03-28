\name{eventsExplorer}
\alias{eventsExplorer}
\title{
  A graphical user interface for defining discrete events in time series
}
\description{
  A graphical user interface for defining discrete events in time
  series. The interface is based on \code{playwith()} from the
  \pkg{playwith} package, which depends on GTK+ being installed on your
  system.
}
\usage{
eventsExplorer(Q, P = NULL,
               ..., type = list(P = "s", "l"), xlab = NULL,
               xscale.components = xscale.components.subticks)
}
\arguments{
  \item{Q, P}{
    \code{\link{zoo}} objects, single time series.
  }
  \item{\dots, type, xlab, xscale.components}{
    further arguments passed to \code{\link{xyplot.ts}}.
  }
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{eventseq}},
  \code{\link{playwith}},
  \code{\link{xyplot.ts}}
}
\examples{
\dontrun{
if (interactive() && require("playwith")) {

data(Queanbeyan)
## wet period
ts74 <- window(Queanbeyan, start = "1974-01-01", end = "1976-12-01")
## dry period
ts93 <- window(Queanbeyan, start = "1993-01-01", end = "1996-01-01")

## To zoom in to a time period: hold Shift and drag horizontally.

eventsExplorer(Q = ts74$Q, P = ts74$P)

#eventsExplorer(ts93$Q)

}}
}
\keyword{ ts }
\keyword{ iplot }
