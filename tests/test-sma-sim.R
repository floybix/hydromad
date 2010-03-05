library(testthat)
library(hydromad)

context("Soil Moisture Accounting models")

set.seed(0)
warmup <- 100
P <- ts(pmax(rnorm(200), 0))
E <- ts(20 + 10 * cos((1:200)/20))
DATA <- cbind(P = P, E = E)

test_that("cmd simulation is the same in R and C", {
    set.seed(0)
    mod0 <- hydromad(DATA, sma = "cmd")
    for (i in 1:5) {
        mod <- simulate(mod0, 1)[[1]]
        hydromad.options(pure.R.code = TRUE)
        pureRsim <- predict(mod)
        hydromad.options(pure.R.code = FALSE)
        expect_that(all(pureRsim >= 0), is_true())
        expect_that(predict(mod), equals(pureRsim))
    }
})

test_that("cwi simulation is the same in R and C", {
    set.seed(0)
    mod0 <- hydromad(DATA, sma = "cwi") 
    for (i in 1:5) {
        mod <- simulate(mod0, 1)[[1]]
        hydromad.options(pure.R.code = TRUE)
        pureRsim <- predict(mod)
        hydromad.options(pure.R.code = FALSE)
        expect_that(all(pureRsim >= 0), is_true())
        expect_that(predict(mod), equals(pureRsim))
    }
})

test_that("bucket simulation is the same in R and C", {
    set.seed(0)
    mod0 <- hydromad(DATA, sma = "bucket") 
    for (i in 1:5) {
        mod <- simulate(mod0, 1)[[1]]
        hydromad.options(pure.R.code = TRUE)
        pureRsim <- predict(mod)
        hydromad.options(pure.R.code = FALSE)
        expect_that(all(pureRsim >= 0), is_true())
        expect_that(predict(mod), equals(pureRsim))
    }
})

test_that("snow simulation is the same in R and C", {
    set.seed(0)
    DATA[,"E"] <- scale(DATA[,"E"]) * 2
    mod0 <- hydromad(DATA, sma = "snow")
    for (i in 1:5) {
        mod <- simulate(mod0, 1)[[1]]
        hydromad.options(pure.R.code = TRUE)
        pureRsim <- predict(mod)
        hydromad.options(pure.R.code = FALSE)
        expect_that(all(pureRsim >= 0), is_true())
        expect_that(predict(mod), equals(pureRsim))
    }
})

test_that("sacramento simulation runs", {
    set.seed(0)
    mod0 <- hydromad(DATA, sma = "sacramento")
    for (i in 1:5) {
        mod <- simulate(mod0, 1)[[1]]
        expect_that(all(predict(mod) >= 0), is_true())
    }
})

test_that("runoffratio simulation runs", {
    set.seed(0)
    mod0 <- hydromad(DATA, sma = "runoffratio")
    for (i in 1:5) {
        mod <- simulate(mod0, 1)[[1]]
        expect_that(all(predict(mod) >= 0), is_true())
    }
})

test_that("dbm simulation runs", {
    set.seed(0)
    mod0 <- hydromad(DATA, sma = "dbm")
    for (i in 1:5) {
        mod <- simulate(mod0, 1)[[1]]
        expect_that(all(predict(mod) >= 0), is_true())
    }
})

test_that("scalar simulation runs", {
    set.seed(0)
    mod0 <- hydromad(DATA, sma = "scalar")
    for (i in 1:5) {
        mod <- simulate(mod0, 1)[[1]]
        expect_that(all(predict(mod) >= 0), is_true())
    }
})
