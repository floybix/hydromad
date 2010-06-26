\name{hydromad.stats}
\alias{hydromad.stats}
\alias{hmadstat}
\title{ Statistics for use in hydromad }
\description{
  These functions provide access to a built-in set of statistics, and
  also allow the user to change or add named statistics.
  Its usage is similar to \code{\link{hydromad.options}}.
}
\usage{
hmadstat(name, DATA = NULL)

hydromad.stats(...)
}
\arguments{
  \item{name}{ character giving the name of a statistic. }
  \item{DATA}{
    in normal usage this is not needed.
    The time series data matrix used to construct the current model.
    If given, any special \code{.()} expressions in the named objective
    function will be evaluated with the given \code{DATA} and \code{Q}.
    These chunks are cached so that repeated evaluation of the objective
    function is faster; but it is then only valid for models using the
    same dataset.
  }
%  \item{negate}{
%    if \code{TRUE}, the returned function is modified to negate the
%    result. This can be useful when defining an objective function to be
%    minimised.
%  }
  \item{...}{
    new stats can be defined, or existing ones modified, using one or
    more arguments of the form 'name = value' or by passing a list of
    such tagged values.  Existing values can be retrieved by supplying
    the names (as character strings) of the components as unnamed
    arguments.
  }
}
\details{
  \code{hmadstat} returns the named statistical function, while
  \code{hydromad.stats} gets or sets a list of named functions.
  
  The default set of available statistics can be listed with
  \code{str(hydromad.stats())}, and consists of:

  \describe{
    \item{abs.err}{
      the mean absolute error.
    }
    \item{RMSE}{
      Root Mean Squared Error.
    }
    \item{bias}{
      bias in data units, \eqn{\sum ( X - Q )}
    }
    \item{rel.bias}{
      bias as a fraction of the total observed flow,
      \eqn{\sum ( X - Q ) / \sum Q}
    }
    \item{r.squared}{
      R Squared (Coefficient of Determination / Nash-Sutcliffe Efficiency),
      \eqn{1 - \sum (Q - X)^2 / \sum (Q - \bar{Q})^2}
    }
    \item{r.sq.sqrt}{
      R Squared (coefficient of determination) using square-root
      transformed data (less weight on peak flows),
      \eqn{1 - \frac{\sum |\sqrt{Q} - \sqrt{X}|^2 }{ \sum |\sqrt{Q} - \bar{\sqrt{Q}}|^2 }}
    }
    \item{r.sq.log}{
      R Squared (coefficient of determination) using log
      transformed data (roughly equal weight to all flow levels),
      \eqn{1 - \frac{\sum |\log{(Q+\epsilon)} - \log{(X+\epsilon)}|^2 }
	            {\sum |\log{(Q+\epsilon)} - \bar{\log{(Q+\epsilon)}}|^2 }}.
      Here \eqn{\epsilon} is the 10 percentile (i.e. lowest
      decile) of the non-zero values of Q.
    }
    \item{r.sq.diff}{
      R Squared (coefficient of determination) of differences between
      successive time steps, i.e. rises and falls.
    }
    \item{r.sq.rank}{
      R Squared (coefficient of determination) using rank-transformed
      data. To avoid trivial differences, data values with an absolute
      difference of less than about 0.1\% of the maximum flow are treated as
      equal, and also those with a relative difference of less than
      about 1\% are treated as equal (technically, the data are transformed as
      \code{round(x, digits = round(3 - log10(max(Q))))} and then
      \code{round(log10(x), digits = 2)} Ties are represented by their
      average rank.
    }
    \item{r.sq.monthly}{
      R Squared (coefficient of determination) with data aggregated into
      calendar months.
    }
    \item{r.sq.smooth7}{
      R Squared (coefficient of determination) using data smoothed with
      a triangular kernel with an effective bandwidth of 7 time steps
      (the actual width of the kernel is about 40\% larger).
    }
    \item{r.sq.seasonal}{
      R Squared (coefficient of determination) where the reference model
      is the mean in each calendar month, rather than the default which
      is the overall mean.
    }
    \item{persistence}{
      R Squared (coefficient of determination) where the reference model
      predicts each time step as the previous observed value. This
      statistic therefore represents a model's performance compared to a
      \emph{naive} one-time-step forecast.
    }
    \item{events.medsums}{
      R Squared (coefficient of determination) of flow totals in and
      between each \emph{event}, defined by observed flow exceeding
      the median level, with a minimum duration of 5 time steps, and
      a minimum separation between events of 5 time steps.  
    }
    \item{events.90sums}{
      R Squared (coefficient of determination) of flow totals in and
      between each \emph{event}, defined by observed flow exceeding
      the 90 percentile level, with a minimum separation between events
      of 5 time steps.  
    }
    \item{events.90max}{
      R Squared (coefficient of determination) of flow maxima in each
      \emph{event}, defined by observed flow exceeding the 90
      percentile level, with a minimum separation between events of 5
      time steps.
    }
    \item{events.90min}{
      R Squared (coefficient of determination) of flow minima
      \strong{between} each \emph{event}, defined by observed flow exceeding
      the 90 percentile level, and with a minimum separation between
      events of 5 time steps. 
    }
    \item{X0}{
      correlation of modelled flow with the model residuals.
    }
    \item{X1}{
      correlation of modelled flow with the model residuals from the
      previous time step. 
    }
    \item{ARPE}{
      Average Relative Parameter Error. Requires that a
      variance-covariance matrix was estimated during calibration.
    }
  }
}
\seealso{
  \code{\link{buildObjectiveFun}},
  \code{\link{fitStat}},
  \code{\link{objFunVal}},
  \code{\link{summary.hydromad}}
}
\examples{
## see current set of stats
names(hydromad.stats())

## extract only one
hmadstat("RMSE")

## calculate stat value with random data
hmadstat("rel.bias")(Q = 1:10, X = 2:11)
}
\keyword{ programming }
