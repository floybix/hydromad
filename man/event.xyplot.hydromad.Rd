\name{event.xyplot.hydromad}
\alias{event.xyplot.hydromad}
\alias{event.xyplot.hydromad.runlist}
\title{ Visualise systematic model errors against covariates }
\description{
  Visualise systematic model errors against covariates.
}
\usage{
## S3 method for class 'hydromad.runlist':
\method{event.xyplot}{hydromad}(x, events,
             formula =
             ~ e(time(Q),first) + e(julian\%\%365,first) +
               sqrt(e(P)) + sqrt(e(U)) +
               sqrt(e(P,max)) + sqrt(e(U,max)) +
               log10(e(Q)+1) + log10(e(lag(Q,-2),first)+1),
             extract = residuals,
             with.U = TRUE,
             ...,
             panel = function(x, y, ..., se) {
                 panel.rug(x = x, alpha = 0.5)
                 panel.abline(h = 0)
                 library(mgcv)
                 panel.smoother(x, y, y ~ s(x), method = "gam", se = se)
                 panel.xyplot(x, y, ...)
             },
             pch = ".", se = TRUE,
             ylab = "residual flow sums in event windows (mm)",
             data = NULL)
}
\arguments{
  \item{x}{
    a \code{hydromad} or \code{hydromad.runlist} object.
  }
  \item{events}{
    event sequence produced by \code{\link{eventseq}}, or a vector
    defining continguous groups on the specified variables.
  }
  \item{formula}{
    formula defining the covariates to plot, as passed to the formula
    method of \code{\link{event.xyplot}}. It may refer to any of the
    variables in the model data frame (\code{observed(x, select = TRUE)});
    additionally it may refer to the value \code{U}, for the effective
    rainfall series derived from the model. 
    Finally \code{\link{julian}} may be referred to (from data time index).
  }
  \item{extract}{
    a function to apply to \code{x} to extract the response variable; by
    default this is \code{residuals} but could be e.g. \code{fitted} or
    \code{function(x) residuals(x, boxcox = TRUE)}.
  }
  \item{with.U}{
    to include modelled effective rainfall \code{U} as a covariate.
  }
  \item{\dots}{
    further arguments passed to \code{\link{event.xyplot}} and on to
    \code{\link{xyplot}}.
  }
  \item{panel, pch, se, ylab}{
    passed to \code{\link{xyplot}} and on to panel functions.
  }
  \item{data}{
    ignored.
  }
}
\value{
  this function returns a trellis object which can be \code{plot}ted.
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{event.xyplot}},
  \code{\link{xyplot}},
  \code{\link{eventapply}}
}
\examples{
data(Cotter)
mod <- hydromad(Cotter[1:1000,], sma = "scalar",
         routing = "armax", rfit = list("sriv", order = c(2,1)))
event.xyplot(mod)
event.xyplot(mod, extract = function(x) residuals(x, boxcox = TRUE))
}
\keyword{ hplot }
\keyword{ ts }

