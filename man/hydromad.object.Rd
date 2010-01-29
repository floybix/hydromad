\name{hydromad.object}
\alias{hydromad.object}
\alias{isValidModel}
\alias{isFullySpecified}
\alias{predict.hydromad}
\alias{print.hydromad}
\alias{summary.hydromad}
\alias{print.summary.hydromad}
\alias{fitted.hydromad}
\alias{residuals.hydromad}
\alias{observed.hydromad}
\title{ Hydromad model object methods }
\description{
  A \code{hydromad} object represents a model, which may be fully
  specified (calibrated) or be defined only by parameter ranges.
  The parameter values are stored along with the observed input and
  output time-series data.
}
\usage{
isValidModel(object)
isFullySpecified(object)

\method{summary}{hydromad}(object, with_uh = TRUE,
	breaks = NULL, coerce = byDays, ...)

\method{predict}{hydromad}(object, newdata = NULL,
	with_uh = TRUE,
	return_state = FALSE,
	return_components = FALSE, ...)
}
\arguments{
  \item{object}{ an object of class \code{hydromad}. }
  \item{breaks}{ if specified, break up the time series and calculate statistics for each chunk.
	Can be a vector of cut points \emph{or} number of intervals to cut into \emph{or} an interval specification,
	one of "sec", "min", "hour", "day", "DSTday", "week", "month" or "year",
  	optionally preceded by an integer and a space, or followed by "s".
  	See \code{\link{cut.Date}}. }
  \item{coerce}{ function to apply to a \code{ts} object to convert its index
	before applying \code{\link{cut}}. Typically \code{\link{byDays}},
  	\code{byYears}, \code{bySecs} or \code{force} (for an integer index). }
  \item{newdata}{ a \code{\link{ts}}-like object or list containing the new time series data.
	It should have named components \code{P}, \code{E} and optionally \code{Q}. }
  \item{with_uh}{ if \code{TRUE}, run the unit hydrograph sub-model;
  	otherwise just return the output from the effective rainfall model. }
  \item{return_state}{ if \code{TRUE}, return state variables from the effective rainfall model. }
  \item{return_components}{ if \code{TRUE}, return a list with each of the modelled flow components
  	(quick flow and slow flow). This is ignored if \code{as.object = TRUE}. }
  \item{\dots}{ further arguments generally passed on to the methods for \code{\link{tf}}. }
}
\details{
  Several standard methods are available for \code{hydromad} objects:
  \code{\link{coef}},
  \code{\link{residuals}},
  \code{\link{fitted}},
  \code{\link{observed}},
  \code{\link{summary}},
  \code{\link{predict}},
  \code{\link{print}},
  etc.

  Plot methods include
  \code{\link{xyplot.hydromad}},
  \code{\link{qqmath.hydromad}},
  \code{\link{errormasscurve}}.

  In any \code{hydromad} object, the time series of effective
  rainfall is available as an item \code{"U"} in the model
  object.
}
\value{
  \code{summary} returns a data frame.

  \code{predict} returns a time series.

  \code{update} returns a new \code{hydromad} object.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{cwi.sim}},
  \code{\link{cmd.sim}},
  \code{\link{bucket.sim}}
}
\examples{

}
\keyword{ methods }

