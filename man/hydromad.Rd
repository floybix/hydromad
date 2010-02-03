\name{hydromad}
\alias{hydromad}
\alias{ihacres}
\alias{update.hydromad}
\title{
  Specify rainfall - runoff (hydrology) models.
}
\description{
  The \code{hydromad} function can be used to specify models with their
  model equations, data, parameters and settings. It allows a general
  two-component structure, where the Soil Moisture Accounting
  (\code{sma}) component and the Routing (\code{routing}) component can
  be arbitrary functions. A method can be specified for fitting the
  dependent routing component.
}
\usage{
hydromad(DATA = zoo(), ...,
        sma = hydromad.getOption("sma"),
        routing = hydromad.getOption("routing"),
        rfit = NULL, weights = NULL,
        warmup = hydromad.getOption("warmup"))

\method{update}{hydromad}(object, ..., newdata,
       sma, routing, rfit, warmup, weights, and.rescale)
}
\arguments{
  \item{DATA, newdata}{
    a \code{\link{ts}}-like object or list with named components:
	\describe{
	  \item{\code{P}}{ time series of areal rainfall depths, usually in mm. }
	  \item{\code{E}}{ time series of potential evapo-transpiration, or more typically,
		temperature as an indicator of this. Optional for some models. }
	  \item{\code{Q}}{ time series of discharge (streamflow) at the catchment outlet.
		This should usually be in units of mm (averaged over the catchment area).
		Use \code{\link{convertFlow}} to convert it. }
	    }
	    }
  \item{\dots}{
    values or ranges for named parameters. Any parameters not given here
    will be taken from defaults given in \code{hydromad.options(sma)}.
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
  The supplied Soil Moisture Accounting (\code{sma}) models are:
  \describe{
    \item{\code{\link{scalar}}}{
      simple linear model: a constant fraction of rainfall reaches
      the stream.
    }
    \item{\code{\link{cwi}}}{
      IHACRES Catchment Wetness Index (CWI) model.
    }
    \item{\code{\link{cmd}}}{
      IHACRES Catchment Moisture Deficit (CMD) model.
    }
    \item{\code{\link{bucket}}}{
      single-bucket models with saturated/unsaturated zones and
      interception. 
    }
    \item{\code{\link{multibucket}}}{
      multiple-bucket models (UNDER DEVELOPMENT).
    }
    \item{\code{\link{sacramento}}}{
      Sacramento Soil Moisture Accounting model (US NWS).
    }
  }

  The supplied Routing (\code{routing}) models are:
  \describe{
    \item{\code{\link{expuh}}}{
      A Unit Hydrograph transfer function defined as a combination of
      exponentially receding components. Each component is defined by
      its time constant and fractional volume, and if there are multiple
      (up to 3) such components they may be arranged in parallel or series.
    }
    \item{\code{\link{uh}}}{
      A general Unit Hydrograph transfer function, defined by
      ARMA (auto-regressive / moving average) coefficients.
    }
  }
  
  The built-in routing fitting (\code{rfit}) methods are:
  \describe{
    \item{\code{\link[tf.ls.fit]{ls}}}{
      Least Squares fitting (with prefiltering and optional
      weighting) of the transfer function. 
    }
    \item{\code{\link[tf.sriv.fit]{sriv}}}{
      Simple Refined Instrumental Variables (SRIV) fitting of the
      transfer function. 
    }
    \item{\code{\link[tf.inverse.fit]{inverse}}}{
      Fitting a transfer function model to the streamflow data,
      independently of the SMA model.
    }
  }

  The built-in model fitting method (for all parameters or just the
  \code{sma} parameters) are:
  \describe{
    \item{\code{\link{fitBySampling}}}{
      returns best result from sampling in parameter ranges using
      random, latin hypercube sampling, or a uniform grid (all combinations). 
    }
    \item{\code{\link{fitByOptim}}}{
      fits using R's \code{\link{optim}} function - which includes a
      Newton type method and a Simplex method. Has multi-start and
      pre-sampling options.
    }
    \item{\code{\link{fitBySCE}}}{
      fit using the Shuffled Complex Evolution (SCE) algorithm.
    }
    \item{\code{\link{fitByDE}}}{
      fit using the Differential Evolution (DE) algorithm.
    }
  }
  
  Methods for working with the model objects are listed under
  \code{\link{hydromad.object}}. 
}
\value{
  the result is a \link[=hydromad.object]{hydromad object}.
}
\references{
  ... TODO
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad.object}}
}
\examples{
data(Cotter)
x <- Cotter[1:1000]
## an unfitted model, with ranges of possible parameter values
modx <- hydromad(x, sma = "cwi", routing = "expuh",
                 tau_s = c(2,100), v_s = c(0,1))
modx
## now try to fit it
fitx <- fitByOptim(modx)
fitx
summary(fitx)
xyplot(fitx)

data(Canning)
x <- window(Canning, start = "1980-01-01", end = "1981-01-01")
## an unfitted model, with ranges of possible parameter values
mod0 <- hydromad(x, sma = "cwi", l = c(0, 100), 
             routing = "uh", rfit = list("ls", order = c(1,1)))
mod0
## now try to fit it
fit1 <- fitByOptim(mod0)
fit1
summary(fit1)
xyplot(fit1)
}
\keyword{ models }
