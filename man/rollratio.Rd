\name{rollratio}
\alias{rollratio}
\alias{easysmooth}
\title{ Ratio of smoothed time series }
\description{
  Smoothe time series. Take the ratio of two smoothed time series.
}
\usage{
easysmooth(data, width, type = 2, na.action = na.pass, ...,
	kern = kernel("fejer", width \%/\% 2, r = type))

rollratio(data, width, ...)
}
\arguments{
  \item{data}{ Multivariate \code{ts} object or \code{zoo} object.
  	The ratio is taken as: the first column divided by second column. }
  \item{width}{ width of the default kernel in time steps. }
  \item{type}{ type of the default ("fejer") kernel: 1 is a simple moving average,
  	2 is a smooth trigonometric function. Other values are also possible;
 	see \code{\link{kernel}}.}
	\item{\dots}{ arguments passed to \code{\link{filter}}:
	allowable are \code{sides = 1} or \code{circular = TRUE}. }
  \item{kern}{ a \code{\link{kernel}} object or vector of filter coefficients.
  	If given this over-rides the other arguments. }
}
\details{
	After smoothing, the ratio is taken as: the first column of \code{data} divided by second column.
}
\value{
	an object like \code{data}.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{ \code{\link[zoo]{rollmean}}, \code{\link{filter}}, \code{\link{kernapply}} }
\examples{
plot(kernel("fejer", 100, r=1))
plot(kernel("fejer", 100, r=2))
plot(kernel("fejer", 100, r=4))
x <- ts(rpois(500, 20))
plot(x)
lines(easysmooth(x, width=50), col=2, lwd=3)
}
\keyword{ ts }
