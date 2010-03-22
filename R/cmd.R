## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' IHACRES Catchment Moisture Deficit (CMD) model.
#' This is CMD soil moisture accounting (effective rainfall) model.
#' It is a conceptual-type model, where input rainfall is
#' partitioned explicitly into drainage, evapo-transpiration, and
#' changes in catchment moisture.
#'
#' The mass balance step is:
#' \deqn{M[t] = M[t-1] - P[t] + E_T[t] + U[t]}
#' 
#' where \eqn{M} represents catchment moisture deficit (CMD),
#' constrained below by 0 (the nominal fully saturated level).
#' P is catchment areal rainfall, \eqn{E_T} is evapo-transpiration, and
#' U is drainage (effective rainfall). All are, typically, in units of
#' mm per time step.
#' 
#' Rainfall effectiveness (i.e. drainage proportion) is
#' a simple \emph{instantaneous} function of the CMD, with a threshold at \eqn{M = d}:
#' \deqn{\frac{\mathrm{d}U}{\mathrm{d}P} = 1 - \min(1, M/d)}{
#'       dU/dP = 1 - min(1, M/d)}
#' 
#' The actual drainage each time step involves the integral of that relation.
#' 
#' Evapo-transpiration (as a proportion of the potential rate, \code{E[t]})
#' is also a simple function of the CMD, with a threshold at \eqn{M = f d}{M = f * d}:
#' \deqn{E_T[t] = e E[t] \min(1, \exp\left(2\left(1 - \frac{M_f}{fd}\right)\right))}{
#'       E_T[t] = e E[t] \min(1, \exp(2(1 - M_f / (fd))))}
#' 
#' Note that the evapo-transpiration calculation is based on \eqn{M_f},
#' which is the CMD after precipitation and drainage have been accounted for.
#'
#' @name IHACRES.CMD.model
#' @aliases IHACRES.CMD.model
#' @param DATA a \code{\link{ts}}-like object or list with named components:
#'   \describe{
#'     \item{\code{P}}{ time series of areal rainfall depths, usually in mm. }
#'     \item{\code{E}}{ time series of potential evapo-transpiration, or more typically,
#'   	temperature as an indicator of this. }
#'     \item{\code{Q}}{ time series of discharge (streamflow) at the catchment outlet.
#'   	This should usually be in units of mm (averaged over the catchment area).
#'   	Use \code{\link{convertFlow}} to convert it. }
#'   }
#' @param d CMD threshold for producing flow (mm).
#' @param f CMD stress threshold as a proportion of \code{d}.
#' @param e temperature to PET conversion factor.
#' @param M_0 starting CMD value (mm).
#' @param return_state to return state variables as well as the effective rainfall.
#' @return \code{cmd.sim} returns the modelled time series of effective rainfall,
#'  or if \code{return_state = TRUE}, a multi-variate time series with named
#'  columns \code{U} (effective rainfall), \code{CMD} and
#'  \code{ET} (evapo-transpiration \eqn{E_T}).
#'
#' @note Normally compiled C code is used for simulation, but if
#' \code{return_state = TRUE} a slower implementation in R is used.
#'
#' @references
#'  Croke, B.F.W. and A.J. Jakeman (2004),
#'  A Catchment Moisture Deficit module for the IHACRES rainfall-runoff model,
#'  \emph{Environmental Modelling and Software}, 19(1): 1-5.
#'  
#'  Croke, B.F.W. and A.J. Jakeman (2005),
#'  Corrigendum to \dQuote{A Catchment Moisture Deficit module for the IHACRES rainfall-runoff model}
#'  [Environ. Model. Softw. 19 (1) (2004) 1-5],
#'  \emph{Environmental Modelling and Software}, 20(7): 977.
#'
#' @seealso \code{\link{hydromad}(sma = "cmd")} to work with models as
#' objects (recommended). 
#'
#' @examples
#'  data(Canning)
#'  x <- cmd.sim(Canning[1:1000,], d = 200, f = 0.7, e = 0.166,
#'               return_state = TRUE)
#'  xyplot(x)
#'
#' @keywords models
#' @export
cmd.sim <-
    function(DATA,
             d, f, e, M_0 = d/2,
             return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    stopifnot(c("P","E") %in% colnames(DATA))
    ## check values
    stopifnot(d >= 0)
    stopifnot(f >= 0)
    stopifnot(e >= 0)
    ## default initial state
    if (is.na(M_0)) M_0 <- d / 2

    ## f is expressed as a proportion of d
    g <- f * d

    P <- DATA[,"P"]
    E <- DATA[,"E"]
    ## skip over missing values (maintaining the state M)
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    ## TODO: return state from C code
    COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (COMPILED) {
        ans <- .C(sma_cmd,
                as.double(P),
                as.double(E),
                as.integer(NROW(DATA)),
                as.double(d),
                as.double(g),
                as.double(e),
                as.double(M_0),
                U = double(NROW(DATA)),
                M = double(NROW(DATA)),
                ET = double(NROW(DATA)),
                NAOK=FALSE, DUP=FALSE, PACKAGE="hydromad")
        U <- ans$U
        M <- ans$M
        ET <- ans$ET
        ## make it a time series object again
        mostattributes(U) <- attributes(DATA)
        class(U) <- "ts"
        attributes(M) <- attributes(ET) <- attributes(U)
    } else {
        ## implementation in R for cross-checking (slow)
        U <- M <- ET <- P
        M_prev <- M_0
        for (t in seq(1, length(P))) {
            ## rainfall reduces CMD (Mf)
            if (P[t] > 0) {
                if (M_prev < d)
                    Mf <- M_prev * exp(-P[t] / d)
                else if (M_prev < d + P[t])
                    Mf <- d * exp((-P[t] + M_prev - d) / d)
                else
                    Mf <- M_prev - P[t]
            } else {
                Mf <- M_prev
            }
            ## drainage (rainfall not accounted for in -dM)
            U[t] <- max(0, P[t] - M_prev + Mf)
            ## evapo-transpiration
            ET[t] <- e * E[t] * min(1, exp(2 * (1 - Mf / g)))
            ET[t] <- max(0, ET[t])
            ## mass balance
            M[t] <- M_prev - P[t] + U[t] + ET[t]
            M_prev <- M[t] <- max(0, M[t])
            ##M[t] <- M_prev <- max(0, Mf + ET[t])
        }
    }
    ## re-insert missing values
    U[bad] <- NA
    if (return_state) return(ts.union(U=U, CMD=M, ET=ET))
    return(U)
}

cmd.ranges <- function()
    list(d = c(50, 550),
         f = c(0.01, 3),
         e = c(0.01, 1.5))
