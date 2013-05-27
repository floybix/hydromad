\name{paretoCatchments}
\alias{paretoCatchments}
\alias{getDominatedLong}
\title{
  Analysis using Pareto-filtering of model performance across simulation
  periods for multiple catchments
}
\description{
  Functions to evaluate model realisations and structures that perform
  well across simulation periods in multiple catchments
}
\usage{
paretoCatchments(stat, csv.name = NA, objectives = "r.squared")
getDominatedLong(res.dom)
}
\arguments{
  \item{stat}{
    data.frame of results, including columns: \code{Catchment}, performance statistics
}
  \item{csv.name}{
    Optional path and prefix in which to write table indicating whether model is
    dominated across periods ("{csv.name}_isdominated_models_catchments.csv") and to
    write a table of non-dominated models ("%s_nondominated_models_catchments.csv")
}
  \item{objectives}{
    Vector of column names containing performance measures to be used by
    \code{\link{areModelsDominated}} to
    determine whether models are dominated across time periods.
  }
  \item{res.dom}{
    data.frame to be converted to long format, from wide format produced
    by \code{\link{areModelsDominated}}
  }
}
\details{
  For each \code{Catchment} in the data.frame of results, models that
  are dominated across simulation periods are identified using
  \code{areModelsDominated}. These results are summarised in a number of ways:

  \itemize{
  \item The number of dominated model realisations in each catchment indicates
  how demanding the Pareto Filtering was.

  \item The number of realisations where each model structure is non-dominated
  suggests the best model structure for each catchment, in terms of its
  ability to perform well across simulation periods.

  \item Table of non-dominated models, in wide format showing
  performance in each simulation period. These models could be used as
  an ensemble of unrejected models.

  \item Minimum performance of all non-dominated models in
  each period. Unacceptably low performance indicates that the model
  selection criteria or model structure do not operate well in
  particular time periods.

  \item Minimum and range of performance of all non-dominated models
  across each period. A large range of non-dominated performance (RNDP)
  indicates variation in performance across periods. The minimum value
  can be used as a threshold to define a feasible set to allow ensemble
  forecasts. If models at or below the minimum value are rejected,
  non-dominated model realisations would be rejected.
}

  \code{getDominatedLong} is a utility function used by \code{paretoCatchments}
}
\value{
  No value returned, used for its side-effects. Prints results described
  in details. Optionally writes csv files, see \code{csv.name}.
}
\author{
Joseph Guillaume
}
\seealso{
  \code{\link{areModelsDominated}}
}
\examples{
data(YeAl97)
paretoCatchments(YeAl97,objectives="E")
}
\keyword{ models }
