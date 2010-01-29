
library(hydromad)

## CANNING RIVER
data(Canning)
## reproduce result from IHACRES Java version tutorial:
cannCal <- window(Canning, start = "1978-01-01", end = "1982-12-31")
## 1. specifying all parameters, including 'c' scale factor
cannRef <- hydromad(cannCal, sma = "cwi", tw = 162, f = 2, l = 300, t_ref = 0, c = 0.000284,
                    routing = "expuh", tau_s = 4.3, delay = 1, warmup = 200)
summary(cannRef)
## 2. setting 'c' parameter from mass balance; different model
cannRef <- hydromad(cannCal, sma = "cwi", tw = 102, f = 1, l = 300, t_ref = 0,
                    routing = "expuh", tau_s = 3.788, delay = 1, warmup = 200)
summary(cannRef)
## 3. calibration with defined parameter ranges
cannMod <- hydromad(cannCal, sma = "cwi", objective = ~ fitStat(Q, X),
                    fit = FALSE, #list(polish = FALSE),
                    routing = "uh",
                    rfit = list("sriv", warmup = 200),
                    tw = seq(2, 200, by = 20),
                    f = seq(0, 10, by = 1),
                    l = seq(0, 400, by = 50),
                    samples = 100,
                    t_ref = 0)
summary(cannMod)
## breakdown of performance in simulation over whole dataset
summary(update(cannMod, newdata = Canning), breaks = "2 years")

## CMD version
cannMod2 <- hydromad(cannCal, sma = "cmd", objective = ~ fitStat(Q, X),
                     routing = "uh",
                     rfit = list("sriv", warmup = 200, order=c(1,1)))
summary(cannMod2)
