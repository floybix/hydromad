## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

gr7j.sim <-
    function(DATA,
             x1, x5, x6, x7, 
             etmult = 1, S_0 = 0.5, 
             return_state = FALSE, 
             mem = 2.0)
{
    
    stopifnot(c("P","E") %in% colnames(DATA))
    ## check values
    stopifnot(x1 >= 0)
    stopifnot(etmult >= 0)
    stopifnot(S_0 >= 0)
    stopifnot(S_0 <= 1)

    inAttr <- attributes(DATA[,1])
    DATA <- as.ts(DATA)
    P <- DATA[,"P"]
    E <- DATA[,"E"] * etmult

    ## skip over missing values (maintaining the state S)
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    
    ans <- .C('sma_gr7j',
                as.double(P),
                as.double(E),
                as.integer(NROW(DATA)),
                as.double(x1),
                as.double(mem),
                as.double(x5),
                as.double(x6),
                as.double(x7),
                as.double(S_0),
                U = numeric(NROW(DATA)),
                S = numeric(NROW(DATA)),
                ET = numeric(NROW(DATA)),
              NAOK=FALSE, DUP=FALSE, PACKAGE="hydromad")
        
        names(ans) <- c("P", "E", "n", "x1", "mem", "x5", "x6", "x7", "S_0", "U", "S", "ET")
        
        U <- ans$U
        S <- ans$S
        ET <- ans$ET
    
    attributes(U) <- inAttr
    ## re-insert missing values
    U[bad] <- NA
    ans <- U
    if (return_state) {
        attributes(S) <- attributes(ET) <- attributes(U)
        ans <- cbind(U=U, S=S, ET=ET)
    }
    return(ans)
    }



gr7j.ranges <- function()
  list(x1 = c(100, 4000), #based on Grigg and Hughes (2018), https://onlinelibrary.wiley.com/doi/full/10.1002/hyp.13282
       x5 = c(0.5, 3.0),
       x6 = c(0.5, 1.5),
       x7 = c(0.01, 2.0),
       # x8 = c(0.3, 1.0),
       etmult = 1)
  

gr7jrouting.sim <-
  function(U, x2, x3, x4, R_0 = 0, split = 0.9,
           return_components = FALSE,
           epsilon = hydromad.getOption("sim.epsilon"),transformed=FALSE)
  {
    if (transformed){
      x2 <- sinh(x2)
      x3 <- exp(x3)
      x4 <- exp(x4)+0.5
    }
    ## check values
    stopifnot(is.numeric(x2))
    stopifnot(x3 >= 0)
    stopifnot(x4 >= 0.5)
    stopifnot(R_0 >= 0)
    stopifnot(R_0 <= 1)
    
    inAttr <- attributes(U)
    U <- as.ts(U)
    
    n <- ceiling(x4)
    m <- ceiling(x4 * 2)
    
    ## S-curves: cumulative proportion of input with time
    n2 <- floor(x4)
    SH1 <- pmin((1:n / x4) ^ (5/2), 1)
    SH2 <- pmin(c(
      0 + 0.5 * (1:n2 / x4) ^ (5/2),
      1 - 0.5 * (2 - n:m / x4) ^ (5/2)),
      1)
    SH2[1:m / x4 > 2] <- 1
    ## unit hydrographs
    UH1 <- diff(c(0, SH1))
    UH2 <- diff(c(0, SH2))
    
    ## skip over missing values (maintaining the state)
    bad <- is.na(U)
    U[bad] <- 0
    
    filter.pad0 <- function(x, f) {
      y <- x
      y[] <- filter(c(rep(0, length(f)), x),
                    filter = f, sides = 1)[-(1:length(f))]
      y
    }
    
    Q9 <- filter.pad0(split * U, UH1)
    Q1 <- filter.pad0((1-split) * U, UH2)
    
    COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (COMPILED) {
      ans <- .C('routing_gr4j',
                as.double(Q9),
                as.double(Q1),
                as.integer(length(U)),
                as.double(x2),
                as.double(x3),
                as.double(R_0),
                Qr = double(length(U)),
                Qd = double(length(U)),
                R = double(length(U)),
                NAOK=FALSE, DUP=FALSE, PACKAGE="hydromad")
      Qr <- ans$Qr
      Qd <- ans$Qd
      R <- ans$R
    } else {
      ## implementation in R for cross-checking (slow)
      Qd <- Qr <- R <- U
      R_prev <- R_0 * x3
      for (t in seq(1, length(U))) {
        Rt_x3 <- R_prev / x3
        ## groundwater exchange term
        F <- x2 * Rt_x3 ^ (7/2)
        ## reservoir level
        R[t] <- max(0, R_prev + Q9[t] + F)
        ## outflow of reservoir
        Qr[t] <- R[t] * ( 1 - (1 + (R[t]/x3)^4)^(-0.25) )
        R[t] <- R[t] - Qr[t]
        ## other store
        Qd[t] <- max(0, Q1[t] + F)
        R_prev <- R[t]
      }
    }
    
    ## zap simulated values smaller than epsilon
    Qr[Qr < epsilon] <- 0
    Qd[Qd < epsilon] <- 0
    
    attributes(Qr) <- attributes(Qd) <- attributes(R) <- inAttr
    
    ## re-insert missing values
    Qr[bad] <- NA
    Qd[bad] <- NA
    
    if (return_components) {
      return(cbind(Xr = Qr, Xd = Qd, R = R))
    } else {
      return(Qr + Qd)
    }
  }

# gr4jrouting.ranges <- function()
#   list(x2 = c(-5, 3),
#        x3 = c(20, 300),
#        x4 = c(1.1, 2.9))

gr7jrouting.ranges <- function()
  list(x2 = c(-10, 10),
       x3 = c(10, 500),
       x4 = c(0.2, 5))


