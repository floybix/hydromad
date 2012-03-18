\name{deconvolution.uh}
\alias{deconvolution.uh}
\title{
  Unit Hydrograph using deconvolution
}
\description{
  Estimates the unit hydrograph using deconvolution.
}
\usage{
deconvolution.uh(P, Q, FWHM = length(P), do.plot = FALSE)
}
\arguments{
  \item{P}{
    Rainfall or effective rainfall time series
}
  \item{Q}{
    Flow time series
}
  \item{FWHM}{
    Full Width Half Maximum to use for Gaussian apodisation function
}
  \item{do.plot}{
    Whether to plot the unit hydrograph immediately
}
}
\details{
  Estimates the unit hydrograph by fourier transform deconvolution of
  the ratio of the cross correlation and auto-correlation functions.
}
\value{
  Vector of length \code{(length(P)-1)*2} representing deconvolved unit
  hydrograph. Note that entries from \code{length(P)} represent are
  obtained by negative lags.
  Optionally plot the unit hydrograph.
}
\references{
  Barry
}
\author{
  Joseph Guillaume
}
\seealso{
  \code{\link{expuh},\link{lambda},\link{armax},\link{powuh}}
}
\examples{
data(Murrindindi)
h <- deconvolution.uh(Murrindindi$P,Murrindindi$Q,do.plot=T)
head(h)
}
\keyword{ ts}
