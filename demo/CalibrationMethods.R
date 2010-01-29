library(hydromad)

data(SalmonBrook)

obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")

P <- obsdat$P
simU <- hydromad.sim(obsdat, sma = "cwi", tw = 30, f = 0.5, c = 1/1000,
                    routing = NULL)
simQ <- hydromad.sim(simU, sma = NULL,
                    routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3)

simQ2 <- hydromad.sim(obsdat, sma = "cwi", tw = 30, f = 0.5, c = 1/1000,
                     routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3)

all.equal(simQ, simQ2)


modDat <- merge(obsdat[,c("P","E")], Q = byDays(simQ))
specs <-
    list(cwiLS21 = hydromad(modDat, sma = "cwi",
             routing = "expuh", rfit = list("ls", order = c(2, 1))),
         cwiSRIV21 = hydromad(modDat, sma = "cwi",
             routing = "expuh", rfit = list("sriv", order = c(2, 1))),
         cwiExpUH = hydromad(modDat, sma = "cwi",
             routing = "expuh", tau_q = c(0,3), tau_s = c(3,100), v_s = c(0,1)))

## objective function surface over parameters
sampLS21 <- simulate(specs$cwiLS21, 200, sampletype = "random", FUN = objFunVal)
sampLS21 <- cbind(attr(sampLS21, "psets"), value = unlist(sampLS21))
levelplot(value ~ tw * f, sampLS21, panel = panel.levelplot.points)
tileplot(value ~ tw * f, sampLS21, aspect = "fill")

## NOTE: sampling is random, so results can vary!
set.seed(0)

#Q <- simulate(specs$cwiExpUH, nsim=1, sampletype = "random", FUN = predict)

fits <- list()
for (i in 1:3) {
    fits[[i]] <-
        with(specs,
             runlist(latinLS = fitBySampling(cwiLS21, samples = 100),
                     latinSRIV = fitBySampling(cwiSRIV21, samples = 100),
                     bfgsLS = fitByOptim(cwiLS21, samples = 64),
                     bfgsSRIV = fitByOptim(cwiSRIV21, samples = 64),
                     simplexLS = fitByOptim(specs$cwiLS21, samples = 64, method = "Nelder-Mead"),
                     simplexSRIV = fitByOptim(specs$cwiSRIV21, samples = 64, method = "Nelder-Mead"),
                     SCE = fitBySCE(cwiExpUH),
                     DE = fitByDE(cwiExpUH))
             )
}

lapply(fits, summary)

lapply(fits, function(x) cbind(summary(x),
                               t(sapply(x, function(x) c(x["funevals"],
                                                         time = x$timing[[1]])))
                               )
       )

lapply(fits, function(fit)
       sapply(fit, function(x) c(x["funevals"], x$timing["elapsed"]))
       )

## TODO: CMD model
## TODO: different simulation parameters
##       - third order routing
##       -
## TODO: noise model
## TODO: tf.invese.fit


modSpec <- hydromad(modDat, sma = "cwi",
                   routing = "expuh", rfit = list("sriv", order = c(2, 1)))

hydromad.options(trace = TRUE)

bfgsFit <- fitByOptim(modSpec, samples = 64, method = "BFGS")
summary(bfgsFit)
