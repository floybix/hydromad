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
  Karel, K. (1990). Membership-set estimation using random scanning and principal component analysis. Mathematics and Computers in Simulation 32(5-6): 535-543 DOI: http://dx.doi.org/10.1016/0378-4754(90)90009-8.

Also possibly related to Iman, R.L. and W.J. Conover (1982). A distribution-free approach to inducing rank correlation among input variables. Communications in Statistics - Simulation and Computation 11(3): 311-334 DOI: http://dx.doi.org/10.1080/03610918208812265.

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