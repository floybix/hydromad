\name{rollccf}
\Rdversion{1.1}
\alias{rollccf}
\alias{xyplot.rollccf}
\title{
  Rolling cross-correlation at given lags.
}
\description{
  Rolling cross-correlation at given lags.
  Can be useful to show how the relationship between two time series
  changes over time, including out-by-one timing errors.
}
\usage{
rollccf(DATA = data.frame(Q = , P = ),
        width = list(365, 90),
        by = 28,
        lags = base.lag + c(0, 1, -1),
        base.lag = estimateDelay(DATA, rises = rises, plot = FALSE),
        rises = FALSE,
        na.action = na.contiguous,
        na.max.fraction = 1/3)

\method{xyplot}{rollccf}(x, data = NULL, ...,
       with.data = TRUE,
       type = list(c("h", "b")),
       type.data = "l",
       par.settings = simpleTheme(pch = ".", cex = 2),
       layout = c(1, length(x$rolls) + with.data * 2),
       strip = strip.default,
       ylim = c(0, 1), xlab = NULL, as.table = TRUE)
}
\arguments{
  \item{DATA}{
    a named list, data frame, time series or zoo object containing the
    two data series.
  }
  \item{width}{
    a list or number specifying the width of window(s), in time steps, in
    which to calculate cross correlations.
  }
  \item{by}{
    temporal resolution: cross correlation is calculated in windows
    spaced every \code{by} time steps.
  }
  \item{lags, base.lag}{
    \code{lags} for which to calculate the cross correlation. By default
    these are based on the overall maximum cross correlation, \code{base.lag}.
  }
  \item{rises}{
    if \code{TRUE}, compute the cross correlation with \emph{rises in}
    streamflow. In this case the streamflow series must be named \code{"Q"}.
  }
  \item{na.action}{
    function to handle missing data in each window (not the whole
    series). This is only applied when the number of missing values is
    less than \code{na.max.fraction}.

    Could be \code{na.exclude}.
  }
  \item{na.max.fraction}{
    if the proportion of missing values in the moving window exceeds
    this value, the corresponding result will be \code{NA}.
  }
  \item{x, data}{
    \code{x} is an object produced by the \code{rollccf}
    function. \code{data} is ignored.
  }
  \item{with.data}{
    if \code{TRUE}, include the original data series in the plot.
  }
  \item{type, type.data}{
    drawing styles for the cross correlation series and input data
    series. See \code{\link[lattice]{panel.xyplot}}.
  }
  \item{par.settings, layout, strip, ylim, xlab, as.table, ...}{
    passed to \code{\link[lattice]{xyplot}}.
  }
}
\details{
  This is a fairly straightforward application of
  \code{\link[zoo]{rollapply}} with the \code{\link{ccf}} function. It
  may be better to do a time-varying regression between the two series.
}
\value{
  \code{rollccf} returns a list of class \code{"rollccf"}, with
  components:
  \item{rolls}{a list of time series of cross correlations.
    One element for each value of \code{width}. }
  \item{data}{time series input data. }
  \item{lags, width, call}{ values of arguments used. }
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{ccf}},
  \code{\link[zoo]{rollapply}}
}
\examples{
data(Canning)
foo <- rollccf(Canning)
xyplot(foo)
}
\keyword{ ts }
