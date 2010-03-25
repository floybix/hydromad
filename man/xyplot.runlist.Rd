\name{xyplot.runlist}
\alias{xyplot.runlist}
\alias{xyplot.runlist}
\alias{qqmath.runlist}
\alias{errormasscurve.runlist}
\title{ Plot results from a set of model runs }
\description{
  Plot results from a set of model runs using Lattice graphics.
}
\usage{
\method{xyplot}{runlist}(x, data,
             residuals. = FALSE,
             coerce = byDays, 
             superpose = FALSE,
             with.P = FALSE, ...,
             x.same = TRUE, y.same = TRUE,
             layout = c(1, NA))

\method{qqmath}{runlist}(x, data,
             residuals=FALSE,
             coerce=byDays, trans=NULL,
             superpose=TRUE,
             auto.key=TRUE,
             type = "b",
             pch=".", cex=2,
             ...)

\method{errormasscurve}{runlist}(x,
	coerce = byDays,
        superpose = TRUE,
        type=c("total", "average", "relative"),
        ...)
}
\arguments{
  \item{x}{
    a \code{runlist} object, which is a list of model objects.
  }
  \item{with.P}{
    to include the input rainfall series in the plot.
  }
  \item{\dots}{
    further arguments are passed on to plotting functions.
  }
}
\details{

}
\value{

}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{runlist}},
  \code{\link[lattice]{xyplot.ts}},
  \code{\link[lattice]{xyplot}},
  \code{\link[lattice]{qqmath}}
}
\examples{
data(HydroTestData)
mod1 <- hydromad(HydroTestData, routing = "expuh",
                 tau_s = 10)
mod2 <- update(mod1, tau_s = 20, tau_q = 5, v_s = 0.5)
mod3 <- update(mod2, lambda = -0.5)
mods <-
  runlist(`single store` = mod1,
          `two stores` = mod2,
          lambda = mod3)
xyplot(mods, with.P = TRUE)
xyplot(mods, scales = list(y = list(log = TRUE)))
}
\keyword{ utilities }


