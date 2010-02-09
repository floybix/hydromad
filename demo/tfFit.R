
library(hydromad)

U <- c(rep(0,10), 1, rep(0,50), 1, rep(0,50)) ## Input signal
X <- filter(2 * U, 0.8, "recursive")          ## Noise free output

ivCoef <- lsCoef <- invCoef <-
    matrix(NA, nrow = 1000, ncol = 2)

hydromad.options(prefilter = 0.5)

pb <- txtProgressBar(max = 1000, style = 3)
for (i in 1:1000) {
   E <- rnorm(length(U)) * 0.02                  ## White noise
   Y <- X + filter(E, c(1.2, -0.8), "recursive") ## Stochastic output
   est1 <- tf.ls.fit(cbind(U=U,Q=Y), normalise=FALSE,
            order = c(n=1,m=0), delay=0, warmup=0)               ## LS estimation
   est2 <- tf.sriv.fit(cbind(U=U,Q=Y), normalise=FALSE,
            order = c(n=1,m=0), delay=0, warmup=0, epsilon=1e-8) ## SRIV estimation
   lsCoef[i,] <- coef(est1, "tau,v")
   ivCoef[i,] <- coef(est2, "tau,v")
   setTxtProgressBar(pb, i)
}
close(pb)
colMeans(ivCoef)
colMeans(lsCoef)


taus <- make.groups(SRIV = ivCoef[,1], 'Least Squares' = lsCoef[,1])
vs <- make.groups(SRIV = ivCoef[,2], 'Least Squares' = lsCoef[,2])
tauplot <- histogram(~ data, taus, subset = (which == "SRIV"), nint = 25)
vplot <- histogram(~ data, vs, subset = (which == "SRIV"), nint = 25)
print(c('Time Constant' = tauplot, 'Steady State Gain' = vplot))

tauplot <- bwplot(which ~ data, taus)
vplot <- bwplot(which ~ data, vs)
print(c('Time Constant' = tauplot, 'Steady State Gain' = vplot))
