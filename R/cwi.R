## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' IHACRES Catchment Wetness Index (CWI) model.
#' This is CWI soil moisture accounting (effective rainfall) model. 
#' Also known as (a part of) the \dQuote{classic} IHACRES model, but
#' includes the extensions for ephemeral catchments of Ye
#' et. al. (1997). It is a metric-type model, where rainfall
#' effectiveness is proportional to a simple antecedent moisture
#' index, and the output is scaled to enforce mass balance.
#'
#' The effective rainfall at each time step is proportional to
#' rainfall, scaled by a soil moisture index \var{s}:
#' \deqn{U_t = c \cdot s_t \cdot P_t}{
#'       U[t] = c * s[t] * P[t]}
#' 
#' Or, if the parameters \code{l} and \code{p} for ephemeral rivers are used:
#' \deqn{U_t = (c (s_t - l))^p \cdot P_t}{
#'       U[t] = (c * (s[t] - l))^p * P[t]}
#' 
#' The soil moisture index \var{s} is calculated by a filter applied to the rainfall,
#' where the drying rate is defined by a \emph{time constant} \eqn{\tau_{\omega,~t}}{tw[t]}:
#' \deqn{s_t = (1 - 1 / \tau_{\omega,t}) s_{t-1} + P_t}{
#'       s[t] = (1 - 1 / tw[t]) * s[t-1] + P[t]}
#' 
#' If \code{f = 0} then the drying time constant is equal to the value of \code{tw}.
#' Otherwise the drying rate varies over time according to the input data \code{E}:
#' \deqn{\tau_{\omega,t} = \tau_\omega \exp(-0.062 f E_t)}{
#'       tw[t] = tw * \exp(- 0.062 * f * E[t])}
#' 
#' Note that the drying rate and effective rainfall are bounded below by 0,
#' a step omitted in the equations above.
#' 
#' @name IHACRES.CWI.model
#' @aliases IHACRES.CWI.model
#' @param DATA a \code{\link{ts}}-like object or list with named components:
#'   \describe{
#'     \item{\code{P}}{ time series of areal rainfall depths, usually in mm. }
#'     \item{\code{E}}{ time series of potential evapo-transpiration, or more typically,
#'   	temperature as an indicator of this. Can be omitted if \code{f = 0}. }
#'     \item{\code{Q}}{ time series of discharge (streamflow) at the catchment outlet.
#'   	This should usually be in units of mm (averaged over the catchment area).
#'   	Use \code{\link{convertFlow}} to convert it. }
#'   }
#' @param tw drying rate at reference temperature (\eqn{\tau_\omega}{tw}).
#'   This is a \emph{time constant}, the number of time steps to reduce
#'   to a fraction \eqn{1/e \approx 37\%}. See definition below. 
#' @param f temperature dependence of drying rate. See definition below.
#'   The case of \code{f=0} describes an invariant drying rate model,
#'   in which case the input data \code{E} is not required.
#' @param c mass balance term.
#'  If this is \code{NA}, it will be set to ensure that \eqn{\sum{U} =
#'  \sum{Q}}, in which case \code{Q} must be included in
#'  \code{DATA}. It can not be omitted for the simulation function
#'  \code{cwi.sim}. 
#' @param l moisture threshold for producing flow (in units of \var{s}).
#'  This can be used together with \code{p} for ephemeral streams.
#' @param p power on soil moisture (above the threshold \code{l}).
#' @param t_ref reference temperature in units of \code{E}
#'  (by convention, 20 deg. C). This is not a parameter; it simply
#'  transforms \code{tw} (scaling \code{tw} by
#'   \code{exp(0.062 * t_ref * f)}.
#' @param s_0 starting value for soil moisture index \var{s}.
#' @param return_state to return state variables as well as the effective rainfall.
#' @return  \code{cwi.sim} returns the modelled time series of effective rainfall,
#'   or if \code{return_state = TRUE}, a multi-variate time series with named
#'   columns \code{U} (effective rainfall), \code{s} (index of soil moisture, \var{s}) and
#'   \code{w} (the recession rate of \var{s}, i.e.
#'   \eqn{(1 - 1 / \tau_{\omega,t})}{
#'     (1 - 1 / tw[t])}.
#'
#' @references
#'  Jakeman, A. J., and G. M. Hornberger (1993),
#'  How much complexity is warranted in a rainfall-runoff model?,
#'  \emph{Water Resources Research}, 29: 2637-2649.
#'  
#'  Jakeman, A.J., I.G. Littlewood, and P.G. Whitehead (1990),
#'  Computation of the instantaneous unit hydrograph and identifiable component flows
#'  with application to two small upland catchments,
#'  \emph{Journal of Hydrology}, 117: 275-300.
#'  
#'  Ye, W., B.C. Bates, N.R. Viney, M. Sivapalan and A.J. Jakeman (1997),
#'  Performance of conceptual rainfall-runoff models in low-yielding ephemeral catchments,
#'  \emph{Water Resources Research}, 33: 153-16.
#'
#' @seealso \code{\link{hydromad}(sma = "cwi")} to work with models as
#' objects (recommended). 
#'
#' @examples
#' data(Canning)
#' x <- cwi.sim(Canning[1:1000,], tw = 162, f = 2, l = 300, t_ref = 0, c = 0.000284,
#'              return_state = TRUE)
#' xyplot(x)
#' 
#' @keywords models
#' @export
cwi.sim <-
    function(DATA,
             tw, f = 0, c,
             l = 0, p = 1,
             t_ref = hydromad.getOption("cwi")$t_ref,
             s_0 = 0,
             return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    if (NCOL(DATA) > 1) stopifnot("P" %in% colnames(DATA))
    if (f > 0) stopifnot("E" %in% colnames(DATA))
    P <- if (NCOL(DATA) > 1) DATA[,"P"] else DATA
    ## special value c = NA used for initial run for scaling
    if (is.na(c))
        c <- 1
    ## check values
    stopifnot(tw >= 0)
    stopifnot(f >= 0)
    stopifnot(c >= 0)
    stopifnot(l >= 0)
    stopifnot(p > 0)
    stopifnot(length(t_ref) == 1)
    ## compute soil moisture index s
    if (f == 0) {
        ## this is an invariant drying rate model
        w <- 1 - 1 / max(tw, 1)
        w <- rep(w, length(P))
        ## use filter_tv rather than filter because it maintains
        ## state over NA gaps (filter resets to 0)
        #s <- filter(P, w, method="recursive", init = s_0)
    } else {
        tw_k <- tw * exp(- 0.062 * f * (DATA[,"E"] - t_ref))
        tw_k <- pmax(tw_k, 1)
        w <- 1 - 1 / tw_k
    }
    s <- filter_tv(P, w, init = s_0)
    ## compute effective rainfall U
    U <- (pmax(s - l, 0) ^ p) * P
    ## TODO: test mass balance scaling with p > 1
    U <- (c ^ p) * U
    if (return_state) return(ts.union(U=U, s=s, w=w))
    return(U)
}

cwi.ranges <- function()
    list(tw = c(0, 100),
         f = c(0, 8),
         c = NA,
         l = 0,
         p = 1,
         t_ref = 20)

absorbScale.hydromad.cwi <- function(object, gain)
{
    #if (!isTRUE(object$estimateScale)) {
    c <- coef(object)[["c"]]
    ## we only want to do this when c is NA (special value)
    if (is.null(c) || !is.na(c))
        return(NULL)
    c <- 1
    p <- coef(object)[["p"]]
    c <- c * (gain ^ p)
    c <- max(c, 0)
    object <- update(object, c = c, and.rescale = FALSE)
    ## TODO - may be significatly faster:
    #object$parlist[["c"]] <- c
    #object$call$c <- signif(c, 6)
    #object$U <- (c ^ p) * object$U
    object
}
