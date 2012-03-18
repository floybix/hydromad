\name{gr4j.transformpar}
\alias{gr4j.transformpar}
\title{
	Transform GR4J parameters
}
\description{
	Apply or reverse transformation of GR4J parameters
}
\usage{
gr4j.transformpar(pars, back = F)
}
\arguments{
  \item{pars}{
  named vector or list of parameters, e.g. as provided by \code{\link{coef.hydromad}}.
}
  \item{back}{
  Whether to transform or untransform (reverse) the parameters.
}
}
\value{
	Named list of transformed/untransformed parameters, depending on value of \code{back}.
}
\author{
	Joseph Guillaume
}
\seealso{
	\code{\link{gr4j}}
}
\examples{
gr4j.transformpar(c(hydromad.getOption("gr4j"),hydromad.getOption("gr4jrouting")))
gr4j.transformpar(c(x1=150,x2=2,x3=50,x4=2),back=F)
}
\keyword{ models }
