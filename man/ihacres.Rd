\name{ihacres}
\alias{ihacres}
\alias{update.ihacres}
\title{
  Specify IHACRES-type models.
}
\description{
  The \code{ihacres} function can be used to specify models with their
  model components, data, parameters and settings. It allows a general
  two-component structure, where the Soil Moisture Accounting
  (\code{sma}) component and the Routing component can be arbitrary
  functions. A method can be specified for fitting the dependent routing
  component.
}
\usage{
ihacres(DATA = xts(), ...,
        sma = "cwi", routing = "expuh", rfit = NULL,
        weights = NULL, warmup = ihacres.getOption("warmup"))

\method{update}{ihacres}(object, ..., newdata,
       sma, routing, rfit, warmup, weights, and.rescale)
}
\arguments{
  \item{DATA}{
    a \code{\link{ts}}-like object or list with named components:
	\describe{
	  \item{\code{P}}{ time series of areal rainfall depths, usually in mm. }
	  \item{\code{E}}{ time series of potential evapo-transpiration, or more typically,
		temperature as an indicator of this. Can be omitted if \code{f = 0}. }
	  \item{\code{Q}}{ time series of discharge (streamflow) at the catchment outlet.
		This should usually be in units of mm (averaged over the catchment area).
		Use \code{\link{convertFlow}} to convert it. }
	    }
	    }
  \item{\dots}{
    values or ranges for named parameters. Any parameters not given here
    will be taken from defaults given in \code{ihacres.options(sma)}.
  }
  \item{sma}{
    name of the Soil Moisture Accounting (SMA) component. May be \code{NULL}, in
    which case the input rainfall will be passed directly to
    \code{routing}. If \code{sma} is specified, a simulation function
    \code{\emph{<sma>}.sim} must exist.
  }
  \item{routing}{
    name of the routing component (i.e. the component which takes in
    effective rainfall from \code{sma} and converts it to streamflow).
    May be \code{NULL}, in which case the
    model output is taken as the output from \code{sma} directly.
  }
  \item{rfit}{
    optional specification for fitting the routing component. If a
    character string is given, then a function
    \code{\emph{<routing>.<rfit>}.fit} must exist.
  }
  \item{weights}{
    weights associated with data points, or a function to apply to the
    data to calculate weights. Not currently used!
  }
  \item{warmup}{
    warmup period in number of time steps.
  }
}
\details{
  not yet.
}
\value{
  the result is an \link[ihacres.object]{\code{ihacres} object}.
}
\references{
  ...
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{ihacres.object}},
  \code{\link{cwi.sim}}
}
\examples{
data(Canning)
x <- Canning["1980"]
mod0 <- ihacres(x, f = 1, sma = "cwi")
mod0
}
\keyword{ models }
