## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' Sacramento Soil Moisture Accounting model.
#' Developed by the US National Weather Service.
#'
#' This description of the model is given by Burnash (1995):
#'
#' \dQuote{The moisture accounting system utilized in the Sacramento
#' Catchment Model is a carefully structured representation of the
#' catchment's soil moisture storage system. It is based on using
#' simple approximations of many of those soil moisture processes
#' which have been reported in the hydrologic literature. The authors
#' have organised these approximations in a manner which would allow
#' the determination of many catchment characteristics from carefully
#' selected portions of the catchment's hydrologic record. Inasmuch as
#' many of the catchment characteristics are related to the soil
#' moisture capabilities of the catchment, an intelligent application
#' of the model start with a good understanding of the three basic
#' types of soil moisture which can potentially influence catchment
#' runoff conditions. These soil moisture types are: (1) Hygroscopic
#' Water, (2) Tension Water and (3) Free Water. } 
#'
#' [...]
#'
#' \dQuote{Streamflow as computed by the Sacramento Catchment Model is
#' the result of processing precipiatation through an algorithm
#' representing the uppermost soil mantle identified as the upper zone
#' and a deeper portion of the soil mantle or lower zone. The
#' algorithm computes runoff in five basic forms. These are (1) direct
#' runoff from permanant and temporary impervious areas, (2) surface
#' runoff due to precipitation occurring at a rate faster than
#' percolation and interflow can take place when both upper zone
#' storages are full, (3) interflow resulting from the lateral
#' drainage of a temporary free water storage, (4) supplemental base
#' flow, and (5) primary base flow.} (Burnash, 1995)
#'
#' @aliases sacramento.sim
#' @param DATA time-series-like object with columns P (precipitation,
#'   mm) and (potential evapo-transpiration, mm). 
#' @param uztwm Upper zone tension water maximum capacity (mm).
#' @param uzfwm Upper zone free water maximum capacity (mm).
#' @param uzk Lateral drainage rate of upper zone free water expressed
#'   as a fraction of contents per day. 
#' @param pctim The fraction of the catchment which produces
#'   impervious runoff during low flow conditions. 
#' @param adimp The additional fraction of the catchment which
#'   exhibits impervious characteristics when the catchment's tension
#'   water requirements are met. 
#' @param zperc Maximum percolation (from upper zone free water into
#'   the lower zone) rate coefficient. 
#' @param rexp An exponent determining the rate of change of the
#'   percolation rate with changing lower zone water contents. 
#' @param lztwm Lower zone tension water maximum capacity (mm).
#' @param lzfsm Lower zone supplemental free water maximum capacity (mm).
#' @param lzfpm Lower zone primary free water maximum capacity (mm).
#' @param lzsk Lateral drainage rate of lower zone supplemental free
#'   water expressed as a fraction of contents per day. 
#' @param lzpk Lateral drainage rate of lower zone primary free water
#'   expressed as a fraction of contents per day. 
#' @param pfree Direct percolation fraction from upper to lower zone
#'   free water (the percentage of percolated water which is available
#'   to the lower zone free water aquifers before all lower zone tension
#'   water deficiencies are satisfied). 
#' @param etmult Multiplier applied to \code{DATA$E} to estimate
#'   potential evapotranspiration. 
#' @param dt Length of each time step in days.
#' @param return_state Not currently supported.
#' @return the simulated effective rainfall (\dQuote{total channel
#'   inflow}), a time series of the same length as the input series.
#'
#' @references
#' Burnash, R.J.C (1995). The NWS River Forecast System -- Catchment Modeling.
#' In: Vijay P. Singh (ed.), \emph{Computer models of watershed hydrology.}
#' Revised edition, Highlands Ranch, Colo. : Water Resources Publications, c1995.
#' \url{http://www.wrpllc.com/books/cmwh.html}.
#'
#' @seealso \code{\link{hydromad}(sma = "sacramento")} to work with models as
#' objects (recommended). 
#'
#' @examples
#'   data(Cotter)
#'   sac0 <- hydromad(Cotter[1:500], sma = "sacramento")
#'   sac0
#'   simulate(sac0, 1)
#' @keywords models
#' @export
sacramento.sim <-
    function(DATA,
             uztwm, uzfwm, uzk, pctim, adimp, zperc, rexp,
             lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree,
             etmult = 1, dt = 1,
             return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    stopifnot(c("P","E") %in% colnames(DATA))
    ## check values
    stopifnot(uztwm >= 0)
    stopifnot(uzfwm >= 0)
    stopifnot(uzk >= 0)
    stopifnot(0 <= pctim && pctim <= 1)
    stopifnot(adimp >= 0)
    stopifnot(zperc >= 0)
    stopifnot(lztwm >= 0)
    stopifnot(lzfsm >= 0)
    stopifnot(lzfpm >= 0)
    stopifnot(lzsk >= 0)
    stopifnot(lzpk >= 0)
    stopifnot(pfree >= 0)
    stopifnot(etmult >= 0)
    stopifnot(dt >= 0)

    xpar <-
        c(uztwm, uzfwm, uzk, pctim, adimp, zperc, rexp,
          lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree)

    P <- DATA[,"P"]
    E <- DATA[,"E"]
    ## skip over missing values
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    ## NOTE: there is not ability to return state variables yet
    ## also no pure.R.code implementation
    U <- .C(sma_sac,
            as.double(P),
            as.double(E),
            as.integer(NROW(DATA)),
            as.double(xpar),
            as.double(etmult),
            as.double(dt),
            U = double(NROW(DATA)),
            NAOK = FALSE, DUP = FALSE, PACKAGE="hydromad")$U
    ## make it a time series object again
    mostattributes(U) <- attributes(DATA)
    class(U) <- "ts"
    ## re-insert missing values
    U[bad] <- NA
    return(U)
}

sacramento.ranges <- function()
    list(uztwm = c(1, 150),
         uzfwm = c(1, 150),
         uzk = c(0.1, 0.5),
         pctim = c(0.00001, 0.1),
         adimp = c(0, 0.4),
         zperc = c(1, 250),
         rexp = c(0, 5),
         lztwm = c(1, 500),
         lzfsm = c(1, 1000),
         lzfpm = c(1, 1000),
         lzsk = c(0.01, 0.25),
         lzpk = c(0.0001, 0.25),
         pfree = c(0, 0.6))
