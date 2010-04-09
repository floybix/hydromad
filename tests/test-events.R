library(testthat)
library(hydromad)

context("Event handling functions")

data(SalmonBrook)
dat <- window(SalmonBrook, start = "1990-01-01", end = "1994-01-01")

evp <- eventseq(dat$P, inter = 7)
evq <- eventseq(dat$Q, thresh = 3)
evq.ts <- eventseq(as.ts(dat$Q), thresh = 3)

test_that("eventseq seems to work", {
    expect_that(evp, is_a("zoo"))
    expect_that(evq.ts, is_a("ts"))
    expect_that(time(evp), equals(time(dat)))
    expect_that(time(evq.ts), equals(time(as.ts(dat))))
    expect_that(coredata(evq), equals(coredata(evq.ts)))
    expect_that(nlevels(evp), equals(39))
    expect_that(nlevels(evq), equals(14))
    expect_that(sum(is.na(coredata(evp))), equals(524))
})

test_that("eventapply seems to work with single series", {
    ## NOTE need to handle functions returning vectors as well as scalars.
    ## (1) scalar result:
    psums <- eventapply(dat$P, evp)
    expect_that(psums, is_a("zoo"))
    expect_that(time(psums), is_a("Date"))
    expect_that(NCOL(psums), equals(1))
    expect_that(NROW(psums), equals(nlevels(evp)))
    ## factor events are not sync'd (cbinded) with the data series
    ## but here we know that they are already synchronised.
    expect_that(psums, is_identical_to(eventapply(dat$P, factor(evp))))
    ## (2) vector (>1) result:
    p2num <- eventapply(dat$P, evp,
                        FUN = function(x) c(mean = mean(x), sd = sd(x)))
    expect_that(colnames(p2num), equals(c("mean", "sd")))
    expect_that(NROW(p2num), equals(nlevels(evp)))
    ## variable length result, with simplify = FALSE
    pvari <- eventapply(dat$P, evp, FUN = coredata, simplify = FALSE)
    expect_that(pvari, is_a("list"))
    expect_that(length(pvari), equals(nlevels(evp)))
    expect_that(names(pvari), equals(format(unname(time(psums)))))
})

test_that("eventapply seems to work with multiple series", {
    ## (1) scalar result with by.column = FALSE
    durs <- eventapply(dat, evp, FUN = nrow, by.column = FALSE)
    expect_that(durs, is_a("zoo"))
    expect_that(NCOL(durs), equals(1))
    expect_that(NROW(durs), equals(nlevels(evp)))
    ## (2) scalar result with by.column = TRUE (the default)
    sums <- eventapply(dat, evp, FUN = sum)
    expect_that(sums, is_a("zoo"))
    expect_that(NCOL(sums), equals(NCOL(dat)))
    expect_that(colnames(sums), equals(colnames(dat)))
    ## (3) vector result with by.column = FALSE
    ## should be exactly the same in this case!
    sums2 <- eventapply(dat, evp, FUN = colSums, by.column = FALSE)
    expect_that(sums, equals(sums2))
    ## (4) vector result with by.column = TRUE
    each2num <- eventapply(dat, evp,
                           FUN = function(x) c(mean = mean(x), sd = sd(x)))
    expect_that(colnames(each2num),
                equals(c("P.mean", "P.sd", "Q.mean", "Q.sd", "E.mean", "E.sd")))
    expect_that(time(each2num), equals(time(sums)))
})


## TODO: test eventinfo
