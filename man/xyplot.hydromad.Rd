\name{xyplot.hydromad}
\alias{xyplot.hydromad}
\alias{qqmath.hydromad}
\alias{tsdiag.hydromad}
\alias{errormasscurve}
\title{ Plot methods for Hydromad model objects }
\description{
  Plot methods...
}
\usage{
\method{xyplot}{hydromad}(x, data = NULL,
	coerce = byDays, all = FALSE,
        superpose = TRUE, with.P = FALSE,
        screens = list(P = "rainfall", "streamflow"),
	...)

\method{qqmath}{hydromad}(x, data,
	coerce = byDays, trans = NULL,
	type = c("p", "l"),
	auto.key = TRUE,
	subset = complete.cases(obsmod),
	...)

errormasscurve(x, ...)
}
\arguments{
  \item{x}{ an object of class \code{hydromad}. }
  \item{data}{ ignored. }
  \item{coerce}{ a function to apply to the time series of observed and modelled values. }
  \item{trans}{ a function to apply to the time series. }
  \item{superpose}{ to overlay time series in one panel. }
  \item{with.P}{ to include the input rainfall series in the plot. }
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
\seealso{ \code{\link{hydromad.object}}, \code{\link[lattice]{xyplot}} }
\examples{

}
\keyword{ hplot }
\keyword{ ts }

