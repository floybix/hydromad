\name{fitStat}
\alias{fitStat}
\title{ Coefficient of Efficiency (R Squared) and other fit stats }
\description{
  Coefficient of Efficiency (R Squared) and other fit stats.
}
\usage{
fitStat(obs, mod, ref = NULL, ..., p = 2,
        trans = NULL, negatives.ok = FALSE,
        na.action = na.pass)
}
\arguments{
  \item{obs}{
    observed data vector.
  }
  \item{mod}{
    model-predicted data vector corresponding to \code{obs}.
  }
  \item{ref}{
    reference model predictions corresponding to \code{obs}. If
    \code{NULL}, \code{ref} is taken as the mean of \code{obs} after
    applying any transformation (\code{trans}).
  }
  \item{...}{
    ignored.
  }
  \item{p}{
    power to apply to absolute residuals (\code{abs(obs - mod)} and
    \code{abs(obs - ref)}.
  }
  \item{trans}{
    a function to apply to each data series before calculating the fit
    statistic. 
  }
  \item{negatives.ok}{
    if \code{FALSE}, the default case, all values in \code{obs},
    \code{mod} and \code{ref} are constrained to be non-negative;
    i.e. negative values are replaced with zero.
  }
  \item{na.action}{
    a function to apply to the time series, which is expected to fill in
    or remove missing values (note, this is optional).
  }
}
\details{
  The result is, after transformation of variables,

  \deqn{sum(abs(obs-mod)^p) / sum(abs(obs-ref)^p)}
  
  If the arguments \code{obs}, \code{mod} or \code{ref} are not plain
  vectors, \code{fitStat} will attempt to merge them together, so
  that corresponding time steps are compared to each other even if the
  time windows are not equal.
}
\value{
  a single numeric value.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad.stats}},
  \code{\link{buildObjectiveFun}},
  \code{\link{objFunVal}},
  \code{\link{summary.hydromad}}
}
\examples{
## generate some data -- it is autocorrelated so the use of these
## stats is somewhat problematic!
set.seed(0)
U <- ts(pmax(0, rgamma(200, shape=0.1, scale=20) - 5))
## simulate error as multiplicative uniform random
Ue <- U * runif(200, min = 0.5, max = 1.5)
## and resample 10 percent of time steps
ii <- sample(seq_along(U), 20)
Ue[ii] <- rev(U[ii])
## apply recursive filter
Q <- filter(U, 0.7, method = "r")
X <- filter(Ue, 0.75, method = "r")

## convert to 'zoo' objects with Date index
Q <- zoo(Q, as.Date("2000-01-01") + 1:200)
X <- zoo(X, time(Q))

xyplot(merge(Q, X), superpose = TRUE)

fitStat(Q, X)

## use absolute residuals rather than squared residuals
## (less weight on large residuals)
fitStat(Q, X, p = 1)

fitStat(Q, X, trans = sqrt)

fitStat(Q, X, trans = function(x) log(x+1))

## use a different reference model (one-step-ahead forecast)
fitStat(Q, X, ref = lag(Q,-1))

## reference as seasonal averages rather than overall average
fitStat(Q, X, ref = ave(Q, months(time(Q))))

## see how the reference model performs in terms of R Squared
fitStat(Q, ave(Q, months(time(Q))))
}
\keyword{ ts }

