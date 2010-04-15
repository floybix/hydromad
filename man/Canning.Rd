\name{Canning}
\alias{Canning}
\docType{data}
\title{
  Rainfall, streamflow and potential evaporation data for Canning River
  at Scenic Drive.
}
\description{
  Daily rainfall, streamflow and potential evaporation for Canning River
  at Scenic Drive upstream of Canning Dam (Western Australia), from 1977-01-01 to 1987-12-31.
  The catchment area is 517 square kilometers.

  The Canning River is a tributary of the Swan River, and is located
  South East of Perth, Western Australia.
}
\usage{data(Canning)}
\format{
  A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.
  It is a regular time series indexed by days, in \code{Date} format.

  There are three columns,
  \code{P} (areal rainfall, mm/day) and
  \code{Q} (streamflow, mm/day).
  \code{E} (potential evapotranspiration, mm/day).
}
\details{
  \describe{
    \item{Rainfall (P)}{
      Daily areal rainfall (mm/day).
      Origin unknown, probably spatial interpolation by Barry Croke.
    }
    \item{Streamflow (Q)}{
      Daily mean streamflow (mm/day).
      Stream gauge ID 616024 "Canning River @ Scenic Drive".
      Latitude -33.4176;
      Longitude 115.9817.
    }
    \item{Potential Evapotranspiration (E)}{
      Origin Unknown.
    }
  }
}
\source{
  Talk to Barry Croke...
}
\examples{
data(Canning)
summary(Canning)
xyplot(Canning)
}
\keyword{datasets}
