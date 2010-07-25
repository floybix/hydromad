\name{xyplot.hydromad}
\alias{xyplot.hydromad}
\alias{xyplot.hydromad.runlist}
\alias{qqmath.hydromad}
\alias{errormasscurve}
\alias{tsdiag.hydromad}
\title{ Plot methods for Hydromad model objects }
\description{
  Plot methods...
}
\usage{
\method{xyplot}{hydromad}(x, data = NULL, ..., scales = list(),
        feasible.bounds = FALSE, col.bounds = "grey80", alpha.bounds = 1,
        border = NA, all = FALSE, superpose = TRUE,
        with.P = FALSE, type = "l",
        type.P = c("h", if ("g" \%in\% type) "g"),
        layout = c(1, NA))

\method{xyplot}{hydromad.runlist}(x, data = NULL, ..., scales = list(),
       all = FALSE, superpose = FALSE,
       with.P = FALSE, type = "l",
       type.P = c("h", if ("g" \%in\% type) "g"),
       layout = c(1, NA))

\method{qqmath}{hydromad}(x, data, ..., all = FALSE, type = "l",
       auto.key = list(lines = TRUE, points = FALSE),
       f.value = ppoints(100), tails.n = 100)

errormasscurve(x, ...)

\method{tsdiag}{hydromad}(object, gof.lag, ...)
}
\arguments{
  \item{x}{ an object of class \code{hydromad}. }
  \item{data}{ ignored. }
  \item{\dots, scales, type, auto.key, layout}{
    further arguments passed on to
    \code{\link{xyplot.zoo}} or \code{\link{qqmath}}.
  } 
  \item{all}{
    passed to \code{fitted()} and \code{observed()}.
  }
  \item{with.P}{
    to include the input rainfall series in the plot.
  }
  \item{type.P}{
    plot type for rainfall, passed to \code{\link{panel.xyplot}}.
  }
  \item{superpose}{
    to overlay observed and modelled time series in one panel.
  }
  \item{f.value, tails.n}{
    arguments to \code{\link{panel.qqmath}}.
  }
  \item{object, gof.lag}{
    passed to the \code{arima} method of \code{\link{tsdiag}}.
  }
}
\value{
  the trellis functions return a trellis object.
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{hydromad.object}},
  \code{\link{xyplot}},
  \code{\link{xyplot.ts}},
  \code{\link{xyplot.list}}
}
\examples{
data(Canning)
cannCal <- window(Canning, start = "1978-01-01", end = "1982-12-31")
mod <-
    hydromad(cannCal, sma = "cwi", tw = 162, f = 2, l = 300,
             t_ref = 0, scale = 0.000284,
             routing = "expuh", tau_s = 4.3, delay = 1, warmup = 200)

xyplot(mod, with.P = TRUE)

c(streamflow = xyplot(mod),
  residuals = xyplot(residuals(mod, type = "h")),
  layout = c(1,2), y.same = TRUE)

xyplot(residuals(mod)) +
  layer(panel.tskernel(..., width = 90, c = 2, col = 1)) +
  layer(panel.tskernel(..., width = 180, c = 2, col = 1, lwd = 2)) +
  layer(panel.tskernel(..., width = 360, c = 2, lwd = 2))

qqmath(mod, scales = list(y = list(log = TRUE)), distribution = qnorm,
       type = c("g", "l"))

errormasscurve(mod, type = c("g", "l"))
}
\keyword{ hplot }
\keyword{ ts }

