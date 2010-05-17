\name{summary.hydromad}
\alias{summary.hydromad}
\alias{print.summary.hydromad}
\alias{print.summaryWithBreaks.hydromad}
\alias{summary.hydromad.runlist}
\title{ Assess and summarise performance of Hydromad models }
\description{
  Assess and summarise performance of Hydromad models.
}
\usage{
\method{summary}{hydromad}(object, breaks = NULL, 
        stats = hydromad.getOption("summary.stats"),
        with.hydrostats = TRUE, na.action = na.exclude, ...)

\method{summary}{hydromad.runlist}(object, ...,
        stats = hydromad.getOption("summary.stats"),
        items = stats)
}
\arguments{
  \item{object}{ an object of class \code{hydromad}. }
  \item{breaks}{
    if specified, break up the time series and calculate statistics for each chunk.
    Can be a vector of cut points \emph{or} number of intervals to cut
    into \emph{or} an interval specification, one of "sec", "min",
    "hour", "day", "DSTday", "week", "month" or "year", optionally
    preceded by an integer and a space, or followed by "s". 
    See \code{\link{cut.Date}}.
  }
  \item{stats}{
    which performance statistics to compute. The definitions are looked
    up by name in \code{hydromad.options("stats")}. See Details.
  }
  \item{with.hydrostats}{
    to also include \code{timesteps}, \code{missing} (number of
    timesteps with missing data), \code{mean.P}, \code{mean.Q},
    \code{runoff} (mean rainfall and streamflow, and their ratio).
  }
  \item{na.action}{
    an optional function to apply to the data to treat missing values before
    calculating statistics.
  }
  \item{\dots}{
    ignored.
  }
  \item{items}{
    passed to \code{\link{summary.runlist}}.
  }
}
\details{
  Definitions of statistics are looked up by name in
  \code{hydromad.options("stats")}, so you can add new statistics
  there. See also \code{\link{objFunVal}}.
  
  The default set of available statistics is:
  \describe{
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
      data. To avoid trivial differences, the data are first rounded to
      \code{(log10(max(Q)) - 3)} decimal places. Ties are represented by
      their minimum rank (so that all zero values will have rank 1).
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
    \item{persistence}{
      R Squared (coefficient of determination) where the reference model
      predicts each time step as the previous observed value. This
      statistic therefore represents a model's performance compared to a
      \dQuote{naive} one-time-step forecast.
    }
    \item{r.sq.seasonal}{
      R Squared (coefficient of determination) where the reference model
      is the mean in each calendar month, rather than the default which
      is the overall mean.
    }
    \item{events.med5sum}{
      R Squared (coefficient of determination) of flow totals in each
      \dQuote{event}, defined by flow exceeding the median observed
      flow, and with a minimum separation between events of 5 time
      steps. Events are defined jointly over the observed and modelled
      time series: events continue while either series exceeds the
      threshold. 
    }
    \item{events.med5max}{
      R Squared (coefficient of determination) of flow maxima in each
      \dQuote{event}, defined by flow exceeding the median observed
      flow, and with a minimum separation between events of 5 time steps.
    }
    \item{events.med5min}{
      R Squared (coefficient of determination) of flow minima
      \strong{between} each \dQuote{event}, defined by flow exceeding
      the median observed flow, and with a minimum separation between
      events of 5 time steps. 
    }
    \item{abs.err}{

    }
    \item{RMSE}{

    }
    \item{U1}{

    }
    \item{X1}{
      correlation of model with error lagged one time step
      \code{cor(head(X,-1), tail(X-Q,-1))}
    }
    \item{ARPE}{
    }
    \item{YIC}{
      
    }
  }
}
\value{
  \code{summary} returns a list with named entries for each of the
  chosen stats (\code{stats}). When \code{breaks} is given it retuns a
  \code{zoo} object of class \code{"summaryWithBreaks"}.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}}
}
\examples{

data(HydroTestData)
mod0 <- hydromad(HydroTestData, sma = "scalar",
                 routing = "expuh", tau_s = 10)
mod0

summary(mod0)

summary(mod0, breaks = "months")
#summary(mod0, breaks = 3)

allstats <- names(hydromad.getOption("stats"))
summary(mod0, stats = allstats)

objFunVal(mod0, objective = hydromad.getOption("stats")$r.sq.log)

str(objFunVal(mod0, objective = hydromad.getOption("stats")))

}
\keyword{ methods }

