\name{xyplot.hydromad}
\alias{xyplot.hydromad}
\alias{qqmath.hydromad}
\alias{errormasscurve}
\alias{tsdiag.hydromad}
\title{ Plot methods for Hydromad model objects }
\description{
  Plot methods...
}
\usage{
\method{xyplot}{hydromad}(x, data = NULL, ..., all = FALSE,
        with.P = FALSE, screens = list(P = "rainfall", "streamflow"),
        type = list(P = "h", "l"), superpose = TRUE)

\method{qqmath}{hydromad}(x, data, ..., all = FALSE,
       type = "l", auto.key = list(lines = TRUE, points = FALSE),
       f.value = ppoints(100), tails.n = 100)

errormasscurve(x, ...)

tsdiag(object, gof.lag, ...)
}
\arguments{
  \item{x}{ an object of class \code{hydromad}. }
  \item{data}{ ignored. }
  \item{\dots, type, auto.key}{
    further arguments passed on to
    \code{\link[zoo]{xyplot.zoo}} or \code{\link[lattice]{qqmath}}.
  } 
  \item{all}{

  }
  \item{with.P}{
    to include the input rainfall series in the plot.
  }
  \item{screens}{
    panels specification, passed to \code{\link[zoo]{xyplot.zoo}}.
  }
  \item{superpose}{
    to overlay time series in one panel.
  }
  \item{f.value, tails.n}{
    arguments to \code{\link[lattice]{panel.qqmath}}.
  }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  the trellis functions return a trellis object.
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{hydromad.object}},
  \code{\link[lattice]{xyplot}},
  \code{\link[lattice]{xyplot.ts}}
}
\examples{
data(Canning)
cannCal <- window(Canning, start = "1978-01-01", end = "1982-12-31")
mod <-
    hydromad(cannCal, sma = "cwi", tw = 162, f = 2, l = 300,
             t_ref = 0, scale = 0.000284,
             routing = "expuh", tau_s = 4.3, delay = 1, warmup = 200)
\dontshow{ .featured_example <- "xyplot" }
xyplot(mod, with.P = TRUE)
\dontshow{ .featured_example <- "qqmath" }
qqmath(mod, scales = list(y = list(log = TRUE)), distribution = qnorm,
       type = c("g", "l"))
\dontshow{ .featured_example <- "errormasscurve" }
errormasscurve(mod, type = c("g", "l"))
}
\keyword{ hplot }
\keyword{ ts }

