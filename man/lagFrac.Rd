\name{lagFrac}
\alias{lagFrac}

\title{
  Non-integer time delay
}
\description{
  Delays P by a non-integer number of timesteps using a filter
}
\usage{
lagFrac(P, TDopt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
    Input to delay (typically rainfall)
}
  \item{TDopt}{
    Positive number with which to delay \code{P}
}
}
\value{
  vector of same length as P
}
\author{
  Thanks to Charles Perrin for original code.
  Implemented in R by Joseph Guillaume
}
\seealso{
  \code{\link{estimateDelayFrac}}
}
\examples{
L<-0.6
P <- c(2,0,5,1,6,10,0,0,0)
V1<-lagFrac(P,L)
}
\keyword{ts}

