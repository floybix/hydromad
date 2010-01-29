\name{estimateDelay}
\alias{estimateDelay}
\title{ Estimate the dead time between input and output }
\description{
  Use cross-correlation to estimate the delay between an input time series and
  (rises in) the corresponding output time series.
}
\usage{
estimateDelay(DATA = data.frame(U=, Q=), rises = TRUE, rank = FALSE,
	n.estimates = 1, lag.max = hydromad.getOption("max.delay"),
        na.action = na.exclude, plot = TRUE, main = NULL, ...)
}
\arguments{
  \item{DATA}{ a \code{\link{ts}}-like object or list with named components:
	\describe{
	  \item{\code{U}}{ input (forcing) time series. }
	  \item{\code{Q}}{ output (response) time series. }
	}
  }
  \item{rises}{ use only rises in the output to estimate delay. }
  \item{rank}{ use rank correlation between input and output. }
  \item{n.estimates}{ number of estimates of delay to produce. }
  \item{lag.max}{ largest delay (in time steps) to consider. }
  \item{na.action}{ handler for missing values.
	The default removes leading and trailing NAs only.
  	Use \code{na.exclude} to remove all NAs, but the result
  	will not be a valid autocorrelation sequence. }
  \item{plot}{ plot the cross-correlation function. }
  \item{main}{ title for plot. }
  \item{\dots}{ further arguments passed to \code{\link{ccf}} or on to \code{\link{plot.acf}}. }
}
\details{
  The estimated delay is the one maximising the cross-correlation function.
}
\value{
  The estimated delay as an integer number of time steps.
  If \code{n.estimates > 1}, that number of integer delays, ordered by the CCF.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{ \code{\link{ccf}} }
\examples{
set.seed(1)
x <- ts(pmax(0, rgamma(200, shape=0.1, scale=20) - 5))
## simulate error as multiplicative uniform random
y <- x * runif(200, min=0.5, max=1.5)
## and resample 10 percent of time steps
ii <- sample(seq_along(y), 20)
y[ii] <- rev(y[ii])
## apply recursive filter and lag
y <- filter(y, 0.8, method="r")
y <- lag(y, -2) # true delay is 2
plot(ts.union(x,y))
estimateDelay(list(U=x, Q=y))
estimateDelay(list(U=x, Q=y), rises=FALSE)
}
\keyword{ ts }

