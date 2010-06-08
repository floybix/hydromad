library(testthat)
library(hydromad)

context("Statistics and summary() methods")

data(SalmonBrook)
obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")

## a fixed model to evaluate
mod <- hydromad(obsdat, sma = "cwi", tw = 30, f = 0.5,
                routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3)

test_that("basic summary() works", {
    s <- summary(mod, with.hydrostats = TRUE)
    expect_that(s, is_a("summary.hydromad"))
    expect_that(is.finite(s$rel.bias), is_true())
    expect_that(print(s), prints_text("Time steps:"))
    s.breaks <- summary(mod, breaks = "12 months")
    expect_that(s.breaks, is_a("zoo"))
})

test_that("all statistics can be evaluated", {
    ss <- objFunVal(mod, hydromad.stats())
    expect_that(ss, is_a("list"))
    expect_that(all(is.finite(unlist(ss))), is_true())
})

test_that("custom objective functions work", {
    spec <- update(mod, v_s = c(0,1))
    set.seed(0)
    fit1 <- fitByOptim1(spec, hmadstat("r.squared", negate = TRUE))
    expect_that(fit1, is_a("hydromad")) 
    set.seed(0)
    fit2 <- fitByOptim1(spec, function(Q,X,...) {
        - hmadstat("r.sq.log")(Q,X) + hmadstat("rel.bias")(Q,X)
    })
    expect_that(fit2, is_a("hydromad"))
    expect_that(coef(fit1)[["v_s"]] != coef(fit2)[["v_s"]], is_true())
})
