\name{YeAl97}
\alias{YeAl97}
\docType{data}
\title{
  Model performance statistics from Ye et al. 1997
}
\description{
  Performance statistics from calibrations of three model structures in
  three catchments in two data periods, assessed using absolute mean deviation, bias and efficiency.
}
\usage{data(YeAl97)}
\format{
  A data frame with 36 observations on the following 7 variables.
  \describe{
    \item{\code{Catchment}}{a character vector, either Canning, Salmon
      or Stones}
    \item{\code{calib.period}}{a character vector, either "First 5Y" or
      "Second 5Y"}
    \item{\code{sim.period}}{a character vector, either "First 5Y" or
      "Second 5Y"}
    \item{\code{Model.str}}{a factor with levels \code{GSFB} \code{IHACRES} \code{LASCAM}}
    \item{\code{A}}{a numeric vector, absolute mean deviation}
    \item{\code{B}}{a numeric vector, bias}
    \item{\code{E}}{a numeric vector, efficiency}
  }
}
\source{
  Ye, W., B.C. Bates, N.R. Viney, M. Sivapalan and A.J. Jakeman (1997). Performance of conceptual rainfall-runoff models in low-yielding ephemeral catchments. \emph{Water Resour. Res.} 33(1): 153-166 DOI: \url{http://dx.doi.org/10.1029/96wr02840}.
}
\references{
 Used by \code{\link{paretoCatchments}}.
}
\examples{
data(YeAl97)
YeAl97
}
\keyword{datasets}
