\name{rotatedSampling}
\alias{rotatedSampling}
\title{
  Sample within rotated feasible parameter space
}
\description{
  Using given feasible parameter set, rotate using eigenvectors and
  sample within rotated hypercube
}
\usage{
rotatedSampling(X, samples, expand = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    matrix or data.frame of feasible parameters
}
  \item{samples}{
    Number of samples to take
}
  \item{expand}{
    Ratio with which to expand bounds in rotated parameter space, as a
    buffer around the original hypercube
}
  \item{\dots}{
    Arguments to pass to \code{\link{parameterSets}}
}
}
\details{
  Parameters may not be in feasible parameter set
}
\value{
  data.frame of parameter values drawn from rotated hypercube
}
\references{
  Karel Keesman
}
\author{
  Joseph Guillaume
}
\seealso{
\code{\link{parameterSets}}
}
\examples{

 X <- matrix(c(1,4,2,5,1,2,2,3),ncol=2)
 Y <- rotatedSampling(X,1e3)
 plot(X[,1],X[,2],col="red",cex=2)
 points(Y[,1],Y[,2],pch=".")
}