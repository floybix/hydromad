library(testthat)
library(hydromad)

set.seed(0)


context("Simulation")

data(SalmonBrook)
obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")

## joint simulation of SMA and routing
simQ <- hydromad.sim(obsdat, sma = "cwi", tw = 30, f = 0.5, c = 1/1000,
                     routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3)

test_that("cwi+expuh simulation looks reasonable", {
    expect_that(NROW(simQ) == NROW(obsdat), is_true())
    expect_that(all(is.finite(simQ)), is_true())
})


context("Calibration results on simulated data")

modDat <- merge(obsdat[,c("P","E")], Q = byDays(simQ))
spec0 <- hydromad(modDat, sma = "cwi", routing = "uh")

test_that("Routing fitting methods work with exact inputs", {
    rspecs <-
        list(ls = update(spec0, rfit = list("ls", order = c(2, 1))),
             sriv = update(spec0, rfit = list("sriv", order = c(2, 1))),
             inv = update(spec0, rfit = list("inverse", order = c(2, 1)))
             )
    ## give actual SMA parameter values, so we are fitting Q from exact U
    fits <- as.runlist(lapply(rspecs, update, tw = 30, f = 0.5, c = 1/1000))
    expect_that(summary(fits$ls)$r.squared > 0.9999, is_true())
    expect_that(summary(fits$sriv)$r.squared > 0.9999, is_true())
    expect_that(summary(fits$inv)$r.squared > 0.9999, is_true())
})

test_that("SMA joint fitting methods work with exact inputs", {
    spec <- update(spec0, routing = "expuh", tau_q = c(0,3), tau_s = c(3,100), v_s = c(0,1))
    expect_that(summary(fitByOptim(spec))$r.squared > 0.99, is_true())
    expect_that(summary(fitBySCE(spec))$r.squared > 0.9999, is_true())
    expect_that(summary(fitByDE(spec))$r.squared > 0.999, is_true())
})
