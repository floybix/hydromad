\name{xyplot.ihacres}
\alias{xyplot.ihacres}
\alias{qqmath.ihacres}
\alias{tsdiag.ihacres}
\alias{errormasscurve}
\title{ Plot methods for IHACRES models }
\description{
  Plot methods...
}
\usage{
\method{xyplot}{ihacres}(x, data,
	coerce = byDays, trans = NULL,
	...)

\method{qqmath}{ihacres}(x, data,
	coerce = byDays, trans = NULL,
	type = c("p", "l"),
	auto.key = TRUE,
	subset = complete.cases(obsmod),
	...)

errormasscurve(x, ...)
}
\arguments{
  \item{x}{ an object of class \code{ihacres}. }
  \item{data}{ ignored. }
  \item{coerce}{ a function to apply to the time series of observed and modelled values. }
  \item{trans}{ a function to apply to the time series. }
  \item{screens}{ panels specification, passed to \code{\link[zoo]{xyplot.zoo}}. }
  \item{col}{ colours. }
  \item{\dots}{ further arguments passed on to \code{\link[zoo]{xyplot.zoo}} or \code{\link[lattice]{qqmath}}. }
  \item{type}{ colours. }
  \item{auto.key}{ colours. }
  \item{subset}{ colours. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
	the trellis functions return a trellis object.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{ \code{\link{ihacres.object}}, \code{\link[lattice]{xyplot}} }
\examples{

}
\keyword{ hplot }
\keyword{ ts }

