\name{eigen.plot}
\alias{eigen.plot}
\title{
Eigenplot
}
\description{
Plot pair-wise eigen values of quadratic response surface model
}
\usage{
eigen.plot(obj, fixed.axis = TRUE)
eigen.plot.single(e, max.value = NA)
}
\arguments{
  \item{obj}{
    \code{\link{rsm}} object, as produced by \code{\link{evalRSM}}
  }
  \item{fixed.axis}{
    Whether to use the same dimensions for each pair of variables. Allows
    comparison of relative width/height of ellipses
  }
  \item{e}{
    List with elements \code{values} and \code{vectors}. Usually determined
    internally by \code{eigen.plot}.
  }
  \item{max.value}{
    If not \code{NA}, \code{xlim} and \code{ylim} are
    \code{c(-max.value,max.value)}, so that separated plots can be
    compared. Usually determined internally by \code{eigen.plot}.
  }
}
\value{
Plot of ellipses for each pair of variables, with larger width/height of
ellipse indicating lower level of identifiability, and rotation of
ellipse away from the vertical indicating level of interaction.
}
\author{
Dario Mavec, Mun-Ju Shin, Joseph Guillaume
}
\seealso{
\code{\link{evalRSM}} which includes an example of \code{eigenplot}
}
\keyword{models}

