library(testthat)
library(hydromad)

set.seed(0)


context("Simulation")

data(SalmonBrook)
obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")

## joint simulation of SMA and routing
simQ <- hydromad.sim(obsdat, sma = "cwi", tw = 30, f = 0.5, c = 1/1000,
                     routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3)

test_that("cwi+expuh(2,1) simulation looks reasonable", {
    expect_that(NROW(simQ) == NROW(obsdat), is_true())
    expect_that(all(is.finite(simQ)), is_true())
})


context("Transfer function parameters and simulation")

test_that("TF parameter conversions work", {
    n0m0 <- c(v_s = 1.1)
    n1m0 <- c(tau_s = 2)
    n1m1 <- c(tau_s = 2, v_s = 0.9)
    n2m0 <- c(tau_s = 30, tau_q = 2, series = 1)
    n2m1 <- c(tau_s = 30, tau_q = 2, v_s = 0.3)
    n2m2 <- c(tau_s = 30, tau_q = 2, v_s = 0.3, v_3 = 0.1)
    ## third order models.
    ## "series = 0"
    ## three components in parallel
    n3s0 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_s = 0.3, v_3 = 0.1)
    ## "series = 1"
    ## two components in series and one in parallel
    ## (s & 3 are in series; q in parallel)
    ## TODO: v_q should be (1 - v_s * v_3) in this case!
    n3s1 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_s = 0.3, v_3 = 0.1, series = 1)
    
    ## "series = 2"
    ## one component in series with two in parallel
    ## (3 in series; s & q in parallel)
    ## TODO: v_q should be (1 - v_s) in this case!
    n3s2 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_s = 0.3, v_3 = 0.1, series = 2)
    ## "series = 3"
    ## three components in series
    n3s3 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, series = 3)
    
    tfConv <- hydromad:::tfParsConvert
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
    U <- ts(pmax(rnorm(100), 0))
    simWith <- function(pars, ...)
        tf.sim(U, pars, ...)

    expect_that(simWith(tfConv(n0m0, "a,b")), equals(simWith(n0m0)))
    expect_that(simWith(tfConv(n1m0, "a,b")), equals(simWith(n1m0)))
    expect_that(simWith(tfConv(n1m1, "a,b")), equals(simWith(n1m1)))
    expect_that(simWith(tfConv(n2m0, "a,b")), equals(simWith(n2m0)))
    expect_that(simWith(tfConv(n2m1, "a,b")), equals(simWith(n2m1)))
    expect_that(simWith(tfConv(n2m2, "a,b")), equals(simWith(n2m2)))
    expect_that(simWith(tfConv(n3s0, "a,b")), equals(simWith(n3s0)))
    expect_that(simWith(tfConv(n3s1, "a,b")), equals(simWith(n3s1)))
    expect_that(simWith(tfConv(n3s2, "a,b")), equals(simWith(n3s2)))
    expect_that(simWith(tfConv(n3s3, "a,b")), equals(simWith(n3s3)))

    xyplot(cbind(simWith(tfConv(n3s0, "a,b")), simWith(n3s0)))

    a <- tfConv(n3s0, "a,b")[(1:3)]
    b <- tfConv(n3s0, "a,b")[-(1:3)]
    xyplot(filter(filter(U, b, sides = 1), a, method="r"))
    xyplot(tf.sim(U, c(a,b))
})

context("UH simulation and inverse simulation")

