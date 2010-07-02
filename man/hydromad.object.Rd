\name{hydromad.object}
\alias{hydromad.object}
\alias{update.hydromad}
\alias{predict.hydromad}
\alias{fitted.hydromad}
\alias{observed.hydromad}
\alias{residuals.hydromad}
\alias{coef.hydromad}
\alias{vcov.hydromad}
\alias{isValidModel}
\alias{print.hydromad}
\title{ Standard methods for Hydromad model objects }
\description{
  A \code{hydromad} object represents a model, which may be fully
  specified (calibrated) or be defined only by parameter ranges.
  The model specification and parameter values are stored along with the
  observed input and output time-series data.
}
\usage{
\method{update}{hydromad}(object, ..., newdata, newpars,
       sma, routing, rfit, warmup, weights, and.rescale)

\method{predict}{hydromad}(object, newdata = NULL,
	which = c("both", "sma", "routing"), ...,
        return_state = FALSE, return_components = FALSE)

\method{fitted}{hydromad}(object, ..., U = FALSE, all = FALSE)

\method{observed}{hydromad}(object, ..., select = "Q", all = FALSE)

\method{residuals}{hydromad}(object, ..., all = FALSE)

\method{coef}{hydromad}(object, which = c("both", "sma", "routing"),
     ..., warn = TRUE, etc = FALSE)

\method{vcov}{hydromad}(object, ...)
%% deviance?

isValidModel(object, ...)
}
\arguments{
  \item{object}{
    an object of class \code{hydromad}.
  }
  \item{\dots}{
    In the \code{update} method, parameter values or ranges for the SMA
    and/or routing simulation functions can be given, as with the
    \code{hydromad()} function. 
    Further arguments to \code{predict} are ignored.
  }
  \item{newdata}{
    a \code{\link{ts}}-like object containing a new time series
    dataset (replacing the original \code{DATA} argument given to the
    \code{hydromad} function).
  }
  \item{newpars}{
    a named list or vector of parameter values; this is equivalent to
    specifying the same values as named arguments (as in
    \dQuote{\dots}).
  }
  \item{sma, routing, rfit, warmup, weights}{
    same arguments as for the \code{\link{hydromad}} function. The
    \code{update} method allows these to be changed on an existing model
    object.     
  }
  \item{and.rescale}{
    set to \code{FALSE} to suppress any automatic adjustment of
    parameters for mass balance.
  }
  \item{which}{
    selects either the SMA or routing model, or both models (the
    default).
    For the \code{predict} method, note that if \code{which =
      "routing"}, then \code{newdata} is treated as the effective rainfall
    (U). 
  }
  \item{return_state}{
    passed to the SMA simulation function, to return state variables.
  }
  \item{return_components}{
    passed to the routing simulation function, to return flow components.
  }
  \item{all}{
    if \code{TRUE}, return the entire time series for which data
    exists. Otherwise, the warmup period (specified in the original call
    to \code{\link{hydromad}} is stripped off.
  }
  \item{U}{
    to return modelled effective rainfall (the output from SMA) rather
    than streamflow.
  }
  \item{select}{
    data series to extract (from the original \code{DATA} argument). Use
    \code{TRUE} to extract all columns.
  }
  \item{warn}{
    by default, \code{coef} gives a warning if the model parameters are
    not fully specifed (i.e. some parameters have ranges rather than
    specific values), because it returns a \code{list} rather than a
    \code{vector} in this case. Setting \code{warn = FALSE} skips the
    warning. 
  }
  \item{etc}{
    by default, \code{coef} returns only the model \emph{parameters}, which are
    defined as being numeric and not wrapped in \code{\link{I}()}. If
    \code{etc = TRUE} is given, then all arguments for the simulation
    function(s) will be returned, which may include other data types
    like logicals or time series. In this case the return value is
    always a \code{list}. 
  }
}
\details{
  Several standard methods are available for \code{hydromad} objects:
  
  (note: these are links to the generic functions only)
  
  \code{\link{update}},
  \code{\link{predict}},
  \code{\link{fitted}},
  \code{\link{observed}},
  \code{\link{residuals}},
  \code{\link{coef}},
  \code{\link{vcov}},
  etc.

  The \code{\link[=summary.hydromad]{summary}} method is documented on a
  different page.
  
  Plot methods include
  \code{\link{xyplot.hydromad}},
  \code{\link{qqmath.hydromad}},
  \code{\link{errormasscurve}}.
  
  \code{isValidModel()} returns \code{TRUE} only if the supplied
  \code{hydromad} object is fully specified and has a calculated output
  series. 
}
\value{
  
  \code{update} returns a new \code{hydromad} object.
  
  \code{predict}, \code{fitted}, \code{observed} and \code{residuals}
  returns time series. 
  
  \code{coef} returns a named numeric vector, or a named \code{list} if
  one or more parameters are not fully specified.
  
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}},
  \code{\link{summary.hydromad}},
  \code{\link{xyplot.hydromad}},
  \code{\link{runlist}}
}
\examples{

}
\keyword{ methods }

