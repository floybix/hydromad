library(testthat)
library(hydromad)

set.seed(0)


context("Transfer Functions")


test_that("cwi+expuh(2,1) simulation looks reasonable", {
    data(SalmonBrook)
    obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")
    ## joint simulation of SMA and routing
    simQ <- hydromad.sim(obsdat, sma = "cwi", tw = 30, f = 0.5, c = 1/1000,
                         routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3)
    expect_that(NROW(simQ) == NROW(obsdat), is_true())
    expect_that(all(is.finite(simQ)), is_true())
})


## Transfer function specifications of various orders
n0m0 <- c(v_s = 1.1)
n1m0 <- c(tau_s = 2)
n1m1 <- c(tau_s = 2, v_s = 0.9)
n2m0 <- c(tau_s = 30, tau_q = 2, series = 1) ## v_q = 1
n2m1 <- c(tau_s = 30, tau_q = 2, v_s = 0.3)
n2m2 <- c(tau_s = 30, tau_q = 2, v_s = 0.3, v_3 = 0.1)
## third order models.
## "series = 0": three components in parallel
n3s0 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_s = 0.3, v_3 = 0.1) ## v_q = 1 - v_s - v_3
## "series = 1": two components in series and one in parallel
## (s & 3 are in series; q in parallel)
n3s1 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_s = 0.3, v_3 = 0.1, series = 1) ## v_q = 1 - (v_s * v_3)
## "series = 2": one component in series with two in parallel
## (3 in series; s & q in parallel)
n3s2 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_s = 0.3, v_3 = 0.1, series = 2) ## v_q = 1 - v_s
## "series = 3": three components in series
n3s3 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_3 = 1, series = 3) ## v_q = 1

tfConv <- hydromad:::tfParsConvert

test_that("TF parameter conversions are consistent", {

    toAbAndBack <- function(tv)
        tfConv(tfConv(tv, "a,b"), "tau,v")[names(tv)]
    expect_that(toAbAndBack(n0m0), equals(n0m0))
    expect_that(toAbAndBack(n1m0), equals(n1m0))
    expect_that(toAbAndBack(n1m1), equals(n1m1))
    expect_that(toAbAndBack(n2m0), equals(n2m0))
    expect_that(toAbAndBack(n2m1), equals(n2m1))
    expect_that(toAbAndBack(n2m2), equals(n2m2))
    ## third-order models.
    expect_that(toAbAndBack(n3s0), equals(n3s0))
    #expect_that(toAbAndBack(n3s1), equals(n3s1))
    expect_that(toAbAndBack(n3s2)[["series"]], equals(2))
    expect_that(toAbAndBack(n3s3), equals(n3s3))

    set.seed(0)
    warmup <- 100
    U <- ts(pmax(rnorm(200), 0))
    armaxSim <- function(pars, ..., init = 0)
        window(armax.sim(U, pars = pars, ..., init = init),
               start = warmup, end = 198)
    expuhSim <- function(pars, ...)
        window(do.call("expuh.sim", c(list(U), as.list(pars), ...)),
               start = warmup, end = 198)

    ## there are small(ish) differences due to initialisation
    tol <- 0.05
    expect_that(armaxSim(tfConv(n0m0, "a,b")), equals(expuhSim(n0m0), tol = tol))
    expect_that(armaxSim(tfConv(n1m0, "a,b")), equals(expuhSim(n1m0), tol = tol))
    expect_that(armaxSim(tfConv(n1m1, "a,b")), equals(expuhSim(n1m1), tol = tol))
    expect_that(armaxSim(tfConv(n2m0, "a,b")), equals(expuhSim(n2m0), tol = tol))
    expect_that(armaxSim(tfConv(n2m1, "a,b")), equals(expuhSim(n2m1), tol = tol))
    expect_that(armaxSim(tfConv(n2m2, "a,b")), equals(expuhSim(n2m2), tol = tol))
    expect_that(armaxSim(tfConv(n3s0, "a,b")), equals(expuhSim(n3s0), tol = tol))
    expect_that(armaxSim(tfConv(n3s1, "a,b")), equals(expuhSim(n3s1), tol = tol))
    expect_that(armaxSim(tfConv(n3s2, "a,b")), equals(expuhSim(n3s2), tol = tol))
    expect_that(armaxSim(tfConv(n3s3, "a,b")), equals(expuhSim(n3s3), tol = tol))
    #xyplot(cbind(armaxSim(tfConv(n2m2, "a,b")), expuhSim(n2m2)), superpose = TRUE)
})


test_that("TF simulation and inverse simulation are consistent", {
    set.seed(0)
    warmup <- 100
    U <- ts(pmax(rnorm(200), 0))
    P <- U + abs(rnorm(U))
    Uw <- window(U, start = warmup, end = end(U)[1]-2)
    Pw <- window(P, start = warmup, end = end(U)[1]-2)
    expuhSim <- function(pars, ...)
        window(do.call("expuh.sim", c(list(U), as.list(pars), ...)),
               start = warmup, end = end(U)[1]-2)
    expect_that(armax.inverse.sim(expuhSim(n0m0), P = Pw, pars = tfConv(n0m0))[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(expuhSim(n1m0), P = Pw, pars = tfConv(n1m0))[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(expuhSim(n1m1), P = Pw, pars = tfConv(n1m1))[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(expuhSim(n2m0), P = Pw, pars = tfConv(n2m0))[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(expuhSim(n2m1), P = Pw, pars = tfConv(n2m1))[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(expuhSim(n2m2), P = Pw, pars = tfConv(n2m2))[-1], equals(Uw[-1]))
    ## third order models seem to need a longer warmup
    set.seed(0)
    warmup <- 200
    U <- ts(pmax(rnorm(300), 0))
    P <- U + abs(rnorm(U))
    Uw <- window(U, start = warmup, end = end(U)[1]-2)
    Pw <- window(P, start = warmup, end = end(U)[1]-2)
    expect_that(armax.inverse.sim(expuhSim(n3s0), P = Pw, pars = tfConv(n3s0))[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(expuhSim(n3s1), P = Pw, pars = tfConv(n3s1))[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(expuhSim(n3s2), P = Pw, pars = tfConv(n3s2))[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(expuhSim(n3s3), P = Pw, pars = tfConv(n3s3))[-1], equals(Uw[-1]))
    #armax.inverse.sim(expuhSim(n3s0), P = Pw, pars = tfConv(n3s0)) - Uw
    #xyplot(cbind(armax.inverse.sim(expuhSim(n3s0), P = Pw, pars = tfConv(n3s0)), Uw), superpose = TRUE)
})
