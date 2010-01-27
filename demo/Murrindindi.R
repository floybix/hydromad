
library(ihacres)

## MURRINDINDI
data(Murrindindi)
## reproduce result from IHACRES Java version tutorial:
murrCal <- window(Murrindindi, start = "1976-06-08", end = "1980-06-08")
## 2. specifying all parameters, including 'c' scale factor
murrRef <- ihacres.cwi(murrCal, tw = 11, f = 1.2, c = 0.003147,
                       uh = tf(pars = list(tau_s = 86.19, tau_q = 2.28, v_s = 0.697),
                                    delay = 0, warmup = 200))
summary(murrRef)
## 2. calibration
murrMod <- ihacres.cwi(murrCal, objective = ~ fitStat(Q, X),
                        fit = list(polish = FALSE),
                        uh = list("sriv", order = c(2,1), warmup = 200))
summary(murrMod)
## breakdown of performance in simulation over whole dataset
summary(update(murrMod, newdata = Murrindindi), breaks = "years")

## CMD version -- try different UH model structures
cmdTest <- tryModelOrders(ihacres.cmd(murrCal), n = 1:3, m = 0:2, delay = 0:1)
ordersSummary(cmdTest)
