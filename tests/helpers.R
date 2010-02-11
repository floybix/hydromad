library(testthat)
library(hydromad)

context("Helper functions")

test_that("flow conversions work correctly on known cases", {
    expect_that(convertFlow(10, from = "cm", to = "in"), equals(3.93700787))
    expect_that(convertFlow(0:10, from = "mm/day", to = "mm/hour"), equals((0:10)/24))
    expect_that(convertFlow(666, from = "m3/sec", to = "ML/day"), equals(666*24*60*60/1000))
    expect_that(convertFlow(4, from = "mm", to = "ML", area.km2 = 2), equals(8))
    expect_that(convertFlow(1, from = "mm", to = "ML"), throws_error())
    expect_that(convertFlow(1, from = "ML / 15 minutes", to = "GL / year"), equals(4*24*365.25/1000))
})

