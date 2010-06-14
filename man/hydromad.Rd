\name{hydromad}
\alias{hydromad}
\alias{ihacres}
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
}
\arguments{
  \item{DATA}{
    a \code{\link{ts}}-like object with named columns:
    \describe{
      \item{\code{P}}{
	time series of areal rainfall depths, usually in mm.
      }
      \item{\code{E}}{
	(optional) time series of potential evapo-transpiration, or more
	typically, temperature as an indicator of this. Required for
	some models but not others.  
      } 
      \item{\code{Q}}{
	(optional) time series of discharge (streamflow) at the
	catchment outlet. Required for calibration but not simulation.
	It should usually be in units of mm (averaged over the catchment
	area). Use \code{\link{convertFlow}} to convert it.
      }
      \item{etc.}{
	other data columns may also be included, and will be accessible
	via the \code{observed()} method.
      }
    }
  }
  \item{\dots}{
    values or ranges for named parameters. Any parameters not given here
    will be taken from defaults given in \code{hydromad.options(sma)}
    and/or \code{hydromad.options(routing)}. In addition, other
    arbitrary arguments may be given here that will be passed on to the
    simulation function(s) and not treated as parameters. To specify a
    numeric object that is not a parameter (such as a time series
    object), wrap it in \code{\link{I}()}. 
  }
  \item{sma}{
    name of the Soil Moisture Accounting (SMA) component. May be
    \code{NULL}, in which case the input rainfall will be passed
    directly to \code{routing}. If \code{sma} is specified, a
    corresponding simulation function \var{sma}\code{.sim} must exist.
  }
  \item{routing}{
    name of the routing component (i.e. the component which takes in
    effective rainfall from \code{sma} and converts it to streamflow).
    May be \code{NULL}, in which case the model output is taken as the
    output from \code{sma} directly. 
  }
  \item{rfit}{
    optional specification for fitting the routing component. If a
    character string is given, then a corresponing function
    \var{routing}\code{.}\var{rfit}\code{.fit} must exist.
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
  The \code{hydromad()} function allows models to be specified with
  the given component models and parameter specifications. The resulting
  object can later be modified using the
  \code{\link[=update.hydromad]{update}} method using the same syntax.
  
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
%    \item{\code{\link{multibucket}}}{
%      multiple-bucket models (UNDER DEVELOPMENT).
%    }
    \item{\code{\link{awbm}}}{
      Australian Water Balance Model: simple 3 bucket model.
    }
    \item{\code{\link{sacramento}}}{
      Sacramento Soil Moisture Accounting model (from US NWS).
    }
    \item{\code{\link{snow}}}{
      simple degree day factor snow model coupled with IHACRES CMD soil
      moisture model. 
    }
    \item{\code{\link{runoffratio}}}{
      rainfall is scaled by the runoff coefficient estimated in a moving
      window (can not be used for prediction).
    }
    \item{\code{\link{dbm}}}{
      rainfall is scaled by corresponding streamflow values raised to a
      power (can not be used for prediction). 
    }
  }

  The supplied Routing (\code{routing}) models are:
  \describe{
    \item{\code{\link{armax}}}{
      A general Unit Hydrograph transfer function, defined by
      ARMA (auto-regressive / moving average) coefficients.
    }
    \item{\code{\link{expuh}}}{
      A Unit Hydrograph transfer function defined as a system of
      exponentially receding components. Each component is defined by
      its time constant and fractional volume, and if there are multiple
      (up to 3) such components they may be in a parallel and/or
      series configuration. 
    }
    \item{\code{\link{lambda}}}{
      a non-linear variant of two exponential stores in parallel: the
      partitioning of flow between the two stores is allowed to vary,
      depending on the magnitude of (effective) rainfall. Intense
      rainfall events tend to be dominated by quick flow, while light
      rainfall events tend to be dominated by slow flow.
    }
%    \item{\code{\link{powuh}}}{
%      [TODO]
%    }
  }
  
  The built-in routing fitting (\code{rfit}) methods are
  (for \code{armax} and \code{expuh} routing):
  \describe{
    \item{\code{\link[=armax.ls.fit]{ls}}}{
      Least Squares fitting (with prefiltering and optional
      weighting) of the transfer function. 
    }
    \item{\code{\link[=armax.sriv.fit]{sriv}}}{
      Simple Refined Instrumental Variables (SRIV) fitting of the
      transfer function. 
    }
    \item{\code{\link[=armax.inverse.fit]{inverse}}}{
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
    \item{\code{\link{fitByDream}}}{
      fit using the DiffeRential Evolution Adaptive Metropolis (DREAM)
      algorithm. This is a Markov Chain Monte Carlo algorithm which also
      gives estimates of the joint probability distribution of
      parameters according to a likelihood function.
    }
  }
  
  Methods for working with the model objects are listed under
  \code{\link{hydromad.object}}. 
}
\value{
  the result from \code{hydromad()} is a
  \link[=hydromad.object]{hydromad object}.
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

## IHACRES CWI model with exponential unit hydrograph
## an unfitted model, with ranges of possible parameter values
modx <- hydromad(x, sma = "cwi", routing = "expuh",
                 tau_s = c(2,100), v_s = c(0,1))
modx
## now try to fit it (roughly - only 1 iteration of SCE!)
fitx <- fitBySCE(modx, control = list(maxit = 1))
fitx
summary(fitx)
xyplot(fitx, with.P = TRUE, type = c("l", "g"))

data(Canning)
x <- window(Canning, start = "1980-01-01", end = "1982-01-01")
xyplot(x)
## IHACRES CWI model with extra parameter l
## Fixed UH (fit once) by inverse method
## an unfitted model, with ranges of possible parameter values
mod0 <- hydromad(x, sma = "cwi", l = c(0, 100), 
             routing = "armax", rfit = list("inverse", order = c(1,1)))
mod0
## now try to fit the free parameters using a Simplex method
fit1 <- fitByOptim(mod0)
fit1
summary(fit1)
xyplot(fit1)
}
\keyword{ models }
