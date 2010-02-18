library(testthat)
library(hydromad)
source("helper-equals-ts.R")
source("helper-expuh-orders.R")

context("Calibration results on simulated data")

data(SalmonBrook)
obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")

## joint simulation of SMA and routing -- the "true" model for testing
simQ <- hydromad.sim(obsdat, sma = "cwi", tw = 30, f = 0.5, c = 1/1000,
                     routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3)

modDat <- merge(obsdat[,c("P","E")], Q = byDays(simQ))
spec <- hydromad(modDat, sma = "cwi", routing = "armax")
## give actual SMA parameter values, so we are fitting Q from exact U
xspec <- update(spec, tw = 30, f = 0.5, c = 1/1000)


test_that("Routing fitting methods work on (2,1) model with exact inputs", {
    fits <-
        runlist(ls = update(xspec, rfit = list("ls", order = c(2,1))),
                sriv = update(xspec, rfit = list("sriv", order = c(2,1))),
                inv = update(xspec, rfit = list("inverse", order = c(2,1)))
                )
    expect_that(summary(fits$ls)$r.squared > 0.9999, is_true())
    expect_that(summary(fits$sriv)$r.squared > 0.9999, is_true())
    expect_that(summary(fits$inv)$r.squared > 0.9999, is_true())
})

U <- fitted(xspec, U = TRUE, all = TRUE)

Q_n0m0 <- expuh.sim(U, pars = n0m0)
Q_n1m0 <- expuh.sim(U, pars = n1m0)
Q_n1m1 <- expuh.sim(U, pars = n1m1)
Q_n2m0 <- expuh.sim(U, pars = n2m0)
Q_n2m1 <- expuh.sim(U, pars = n2m1)
Q_n2m2 <- expuh.sim(U, pars = n2m2)
Q_n3s0 <- expuh.sim(U, pars = n3s0)
Q_n3s1 <- expuh.sim(U, pars = n3s1)
Q_n3s2 <- expuh.sim(U, pars = n3s2)
Q_n3s3 <- expuh.sim(U, pars = n3s3)

test_that("Least squares (armax) fitting works for all orders with exact inputs", {
    ## note hydromad.getOption("warmup") == 100
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n0m0), order = c(0,0))), equals_in_window(Q_n0m0, tol = 1e-5))
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n1m0), order = c(1,0))), equals_in_window(Q_n1m0, tol = 1e-5))
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n1m1), order = c(1,1))), equals_in_window(Q_n1m1, tol = 1e-5))
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n2m0), order = c(2,0))), equals_in_window(Q_n2m0, tol = 1e-5))
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n2m1), order = c(2,1))), equals_in_window(Q_n2m1, tol = 1e-5))
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n2m2), order = c(2,2))), equals_in_window(Q_n2m2, tol = 1e-5))
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n3s0), order = c(3,2))), equals_in_window(Q_n3s0, tol = 1e-5))
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n3s1), order = c(3,2))), equals_in_window(Q_n3s1, tol = 1e-5))
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n3s2), order = c(3,1))), equals_in_window(Q_n3s2, tol = 1e-5))
    expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n3s3), order = c(3,0))), equals_in_window(Q_n3s3, tol = 1e-5))
})

test_that("SRIV (armax) fitting works for all orders with exact inputs", {
    ## note hydromad.getOption("warmup") == 100
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n0m0), order = c(0,0))), equals_in_window(Q_n0m0, tol = 1e-5))
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n1m0), order = c(1,0))), equals_in_window(Q_n1m0, tol = 1e-5))
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n1m1), order = c(1,1))), equals_in_window(Q_n1m1, tol = 1e-5))
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n2m0), order = c(2,0))), equals_in_window(Q_n2m0, tol = 1e-5))
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n2m1), order = c(2,1))), equals_in_window(Q_n2m1, tol = 1e-5))
    ## TODO: there's an NA at the end of the sriv fitted from here down
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n2m2), order = c(2,2))), equals_in_window(Q_n2m2, tol = 1e-5))
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n3s0), order = c(3,2))), equals_in_window(Q_n3s0, tol = 1e-5))
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n3s1), order = c(3,2))), equals_in_window(Q_n3s1, tol = 1e-5))
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n3s2), order = c(3,1))), equals_in_window(Q_n3s2, tol = 1e-5))
    expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n3s3), order = c(3,0))), equals_in_window(Q_n3s3, tol = 1e-5))
})

test_that("Inverse (armax) fitting works for all orders with exact inputs", {
    mod <- update(xspec, routing = "expuh")
    U <- fitted(xspec, U = TRUE, all = TRUE)
    
    Q <- expuh.sim(U, pars = n0m0)
    fit <- armax.ls.fit(cbind(U = U, Q = Q), order = c(0,0))
    expect_that(fitted(fit), equals_in_window(Q))

    ## TODO
                                        #mod1 <- do.call("update", c(list(mod), as.list(n0m0)))
    mod1 <- mod
    coef(mod1) <- c(coef(xspec), as.list(n0m0))
    sim1 <- fitted(mod1, all = TRUE)
    
    fit <- update(mod, rfit = list("ls", order = c(0,0)))
    stopifnot(summary(fit)$r.squared > 0.999, is_true())

})

test_that("SMA joint fitting methods work with exact inputs", {
    set.seed(0)
    jspec <- update(spec, routing = "expuh", tau_q = c(0,3), tau_s = c(3,100), v_s = c(0,1))
    ## TODO: suppress optim() output
    expect_that(summary(fitByOptim(jspec))$r.squared > 0.97, is_true())
    expect_that(summary(fitBySCE(jspec))$r.squared > 0.9999, is_true())
    expect_that(summary(fitByDE(jspec))$r.squared > 0.999, is_true())
})

