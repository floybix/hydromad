\name{byDays}
\alias{byDays}
\alias{byYears}
\alias{bySecs}
\title{ Convenient time series conversion }
\description{
  These convenience functions are based on \code{\link[zoo]{as.zoo}}.
  They simply coerce the \code{\link[zoo]{index}} to 
  \code{\link{Date}} (\code{byDays}),
  \code{\link[zoo]{yearmon}} (\code{byYears}) or 
  \code{\link{POSIXct}} (\code{bySecs}).
}
\usage{
byDays(x)
byYears(x)
bySecs(x, tz = "")
}
\arguments{
  \item{x}{ an object to be coerced to \code{\link[zoo]{zoo}}, typically a \code{\link{ts}} object. }
  \item{tz}{ time zone for printing. }
}
\value{
  an object of class \code{zoo}, with the index class as described above
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{ \code{\link[zoo]{index}}, \code{\link[zoo]{zoo}} }
\examples{
x <- ts(1:10, start=2000)
time(x)
time(byDays(x))
time(byYears(x))
time(bySecs(x))

## byYears can be used for a monthly index in fractions of years
x <- ts(1:10, start=2000, frequency=4)
time(x)
time(byYears(x))

## bySecs can be used for a POSIXct index
tt <- seq(as.POSIXct("2006-08-02 15:00"), by="hours", len=10)
x <- zoo(1:10, order.by=tt)
(xf <- filter(x, c(0.5,0.5))) # was converted to ts
bySecs(xf) # convert back
}
\keyword{ ts }

