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
    ss <- objFunVal(mod, hydromad.getOption("stats"))
    expect_that(ss, is_a("list"))
    expect_that(all(is.finite(unlist(ss))), is_true())
})
