\name{event.xyplot}
\alias{event.xyplot}
\alias{event.xyplot.formula}
\title{ Scatterplots with variables aggregated in event windows }
\description{
  Scatterplots with variables aggregated in event windows.
}
\usage{
event.xyplot(x, ...)

\method{event.xyplot}{formula}(x, data = list(), events, ...,
             response = NULL, eFUN = sum,
             as.table = TRUE, layout = NULL,
             auto.key = NCOL(response) > 1,
             xlab = "event-aggregated covariate values",
             ylab = "event-aggregated response")
}
\arguments{
  \item{x}{
    an object for the generic method; in the \code{formula} method this
    is a \code{\link{formula}} as described in the Details section.
    All terms in the formula should be aggregated in the given
    \code{events}; for this purpose a special function \code{e()} is
    defined as a shorthand for \code{\link{eventapply}}, with arguments
    \code{e(X, FUN = sum, ...)}.
    The formula may specify a response (the left hand side of formula),
    which is treated as a single term, although it can be a matrix with
    multiple columns. Alternatively, the response matrix can be given
    in \code{response}.
    The covariate terms are also allowed to be matrix-like with the same
    number of columns as the response, in which case each column is
    plotted against the corresponsing column of the response.
  }
  \item{data}{
    data source giving variables used in the formula. Typically a
    \code{data.frame}. Passed to \code{\link{model.frame}}. 
  }
  \item{events}{
    event sequence produced by \code{\link{eventseq}}, or a vector
    defining continguous groups on the specified variables.
  }
  \item{\dots}{
    further arguments passed to \code{\link{xyplot}} (actually
    \code{\link{xyplot.list}}). 
  }
  \item{response}{
    a response vector or matrix corresponding to the variables in the
    formula. If this is given, any response in the formula is ignored.
  }
  \item{eFUN}{
    a function to aggregate \code{response}, if it is given, in
    the given \code{events}. If a reponse appears in the formula it
    should be in a \code{e()} term, and \code{eFUN} is not used.
  }
  \item{as.table, layout, auto.key, xlab, ylab}{
    passed to \code{\link{xyplot}}.
  }
}
\value{
  this function returns a trellis object which can be \code{plot}ted.
}
\author{
  Felix Andrews \email{felix@nfrac.org}
}
\seealso{
  \code{\link{xyplot.list}},
  \code{\link{xyplot}},
  \code{\link{panel.xyplot}},
  \code{\link{eventapply}}
}
\examples{
data(Canning)
ev <- eventseq(Canning$P, thresh = 20, inthresh = 1, indur = 3,
               continue = FALSE)
event.xyplot(log10(e(Q,sum)+0.01) ~ e(P,sum) + sqrt(e(P,max)) + log10(e(lag(Q,-1),first)+0.01),
             data = Canning, events = ev, 
             ylab = "event-aggregated flow sums (mm, log10)")
## change to nicer strip labels and add a loess smoother
foo <- trellis.last.object(type = c("p", "smooth"))
dimnames(foo)[[1]] <-
    paste("vs.", c("total rain (mm)", "max rain (sqrt mm/day)", "log10 antecedent flow (mm/day)"))
foo

## multiple response variables
zquantile <- function(x,...) quantile(coredata(x),...)
event.xyplot(cbind(e(Q,zquantile,0.9),e(Q,zquantile,0.98)) ~
             e(P,max) + e(P,mean) + e(P,sum),
             data = Canning, events = ev,
             scales = list(log = TRUE),
             xscale.components = xscale.components.log10ticks,
             yscale.components = yscale.components.log10.3)
}
\keyword{ hplot }
\keyword{ ts }

