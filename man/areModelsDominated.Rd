\name{areModelsDominated}
\alias{areModelsDominated}
\title{
  Are Models Dominated?
}
\description{
  Identify dominated model realisations, that are inferior to another
  model in all simulation periods.
}
\usage{
areModelsDominated(res, objectives = "r.squared")
}
\arguments{
  \item{res}{
    data.frame of results, including the column \code{sim.period} and
    the columns named in \code{objectives}. At least one of the following column names
    should be included as id variables:
    \code{Model.str},\code{Catchment},\code{calib.period},\code{Cal.objfn}
    Other columns will be ignored.
  }
  \item{objectives}{
    Vector of column names containing performance measures. We assume
    higher values are better. Values should be transformed prior to use.
  }
}
\details{
  The data.frame of results \code{res} is transformed to wide format with a column
  for performance in each \code{sim.period} for each of the
  \code{objectives}.
  Each row is interpreted as a candidate model realisation.
  The resulting matrix is filtered using \code{\link{paretoFilter}} to
  identify dominated models.
}
\value{
  A wide-format data.frame with id variable columns, a column for each
  \code{sim.period} value, and a column \code{dominated} indicating
  whether another model is better in all simulation periods.
}
\author{
Joseph Guillaume
}
\seealso{
  \code{\link{paretoFilter}}, \code{\link{paretoCatchments}}
}
\examples{

data(YeAl97)

## Dataset consisting of results for two simulation periods,
##  obtained by calibration in the same periods with different
##  model structures.
results <- subset(YeAl97,Catchment=="Salmon")

## Table indicating whether models defined by their
##  calib.period and Model.str are dominated according
##  to the objective E
areModelsDominated(results,objectives="E")
}
\keyword{ models }