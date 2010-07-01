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
    If given, the returned statistic function will be pre-evaluated with
    the given data matrix (which must have a column \code{Q} of observed
    data). Technically, any special \code{.()} expressions in the named objective
    function will be evaluated with the given \code{DATA} and \code{Q}.
    These chunks are cached so that repeated evaluation of the objective
    function is faster; but it is then only valid for models using the
    same dataset.
    In normal usage this is not needed.
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
      R Squared (Nash-Sutcliffe Efficiency),
      \eqn{1 - \sum (Q - X)^2 / \sum (Q - \bar{Q})^2}
    }
    \item{r.sq.sqrt}{
      R Squared using square-root
      transformed data (less weight on peak flows),
      \eqn{1 - \frac{\sum |\sqrt{Q} - \sqrt{X}|^2 }{ \sum |\sqrt{Q} - \bar{\sqrt{Q}}|^2 }}
    }
    \item{r.sq.log}{
      R Squared using log
      transformed data (roughly equal weight to all flow levels),
      \eqn{1 - \frac{\sum |\log{(Q+\epsilon)} - \log{(X+\epsilon)}|^2 }
	            {\sum |\log{(Q+\epsilon)} - \bar{\log{(Q+\epsilon)}}|^2 }}.
      Here \eqn{\epsilon} is the 10 percentile (i.e. lowest
      decile) of the non-zero values of Q.
    }
    \item{r.sq.diff}{
      R Squared using differences between
      successive time steps, i.e. rises and falls.
    }
    \item{r.sq.monthly}{
      R Squared with data aggregated into
      calendar months.
    }
    \item{r.sq.smooth7}{
      R Squared using data smoothed with
      a triangular kernel with an effective bandwidth of 7 time steps
      (the actual width of the kernel is about 40\% larger).
    }
    \item{r.sq.seasonal}{
      R Squared where the reference model
      is the mean in each calendar month, rather than the default which
      is the overall mean.
    }
    \item{persistence}{
      R Squared where the reference model
      predicts each time step as the previous observed value. This
      statistic therefore represents a model's performance compared to a
      \emph{naive} one-time-step forecast.
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
  \code{\link{nseStat}},
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
