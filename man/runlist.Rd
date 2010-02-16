\name{runlist}
\alias{runlist}
\alias{as.runlist}
\alias{[.runlist}
\alias{xyplot.runlist}
\alias{qqmath.runlist}
\alias{errormasscurve.runlist}
\alias{print.summary.runlist}
\title{ Work with a set of model runs }
\description{
  A \code{runlist} object is simply a list of model objects.
}
\usage{
runlist(...)

as.runlist(x, ...)

\method{xyplot}{runlist}(x, data,
             residuals. = FALSE,
             coerce = byDays, 
             superpose = FALSE,
             with.P = FALSE, ...,
             x.same = TRUE, y.same = TRUE,
             layout = c(1, NA))

\method{qqmath}{runlist}(x, data,
             residuals=FALSE,
             coerce=byDays, trans=NULL,
             superpose=TRUE,
             auto.key=TRUE,
             type = "b",
             pch=".", cex=2,
             ...)

\method{errormasscurve}{runlist}(x,
	coerce = byDays,
        superpose = TRUE,
        type=c("total", "average", "relative"),
        ...)
}
\arguments{
  \item{\dots}{ for \code{runlist}, a named list of model objects,
    specified as in \code{\link{list}}.
    For plot methods, further arguments are passed on to plotting
    functions. }
  \item{x}{
    for \code{as.runlist}, a simple \code{list}; in other cases, a
    \code{runlist} object. }
  \item{with.P}{ to include the input rainfall series in the plot. }  
}
\details{

}
\value{

}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{ \code{\link[zoo]{xyplot.zoo}}, \code{\link[lattice]{qqmath}} }
\examples{

}
\keyword{ utilities }


