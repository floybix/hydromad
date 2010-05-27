library(testthat)
#library_if_available(hydromad)
library(hydromad)


context("Calibration results with real data")


test_that("cwi+expuh optim gives reasonable result", {
    data(Cotter)
    x <- Cotter[1:1000]
    modx <- hydromad(x, sma = "cwi", routing = "expuh",
                     tau_s = c(1, 100), v_s = c(0, 1))
    ## now try to fit it
    set.seed(0)
    fitx <- fitByOptim(modx)
    s <- summary(fitx)
    expect_that(s$r.squared > 0.65, is_true())
    expect_that(s$r.sq.log > 0.8, is_true())
    expect_that(abs(s$rel.bias) < 0.01, is_true())
})


test_that("cwi+expuh optim/SRIV gives reasonable result", {
    data(Cotter)
    x <- Cotter[1:1000]
    modx <- hydromad(x, sma = "cwi", routing = "expuh",
                     rfit = list("sriv", order = c(2,1)))
    ## now try to fit it
    set.seed(0)
    fitx <- fitByOptim(modx)
    s <- summary(fitx)
    expect_that(s$r.squared > 0.65, is_true())
    expect_that(s$r.sq.log > 0.8, is_true())
    expect_that(abs(s$rel.bias) < 0.01, is_true())
})


test_that("cmd+powuh SCE gives reasonable result", {
    data(Cotter)
    x <- Cotter[1:1000]
    modx <- hydromad(x, sma = "cmd", routing = "powuh")
    ## now try to fit it
    set.seed(0)
    fitx <- fitBySCE(modx)
    s <- summary(fitx)
    expect_that(s$r.squared > 0.7, is_true())
    expect_that(s$r.sq.log > 0.8, is_true())
    expect_that(abs(s$rel.bias) < 0.1, is_true())
})

