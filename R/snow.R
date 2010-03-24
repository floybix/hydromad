## hydromad: Hydrological Modelling and Analysis of Data


# Simple degree day factor snow model with IHACRES CMD soil moisture model
## coded by Jarkko Koskela @tkk.fi 2010-02-26
snow.sim <-
    function(DATA, Tmax, Tmin, kd, kf, rcap, Tmelt = Tmin,
             cr = 1, cs = 1, LSWE_0 = 0, ISWE_0 = 0,
             d, f, e, M_0 = d/2, return_state = FALSE)
{
    ## get data into the right form
    DATA <- as.ts(DATA)
    stopifnot(c("P","E") %in% colnames(DATA))
    ## check values
    stopifnot(0 <= kd)
    stopifnot(0 <= kf)
    stopifnot(0 <= cr)
    stopifnot(0 <= cs)
    stopifnot(0 <= rcap)
    Tmin <- min(Tmax, Tmin)
    P <- DATA[,"P"]
    E <- DATA[,"E"]

    ## rainfall or snowfall
    fr <- (E - Tmin) / (Tmax - Tmin)
    fr <- pmax(pmin(fr, 1), 0)
    Prain <- fr * cr * P
    Psnow <- (1-fr) * cs * P

    COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (COMPILED) {
        ans <- .C(sma_snow,
                  as.double(Prain),
                  as.double(Psnow),
                  as.double(E),
                  as.integer(NROW(DATA)),
                  as.double(kd),
                  as.double(kf),
                  as.double(rcap),
                  as.double(Tmelt),
                  as.double(LSWE_0),
                  as.double(ISWE_0),
                  U = double(NROW(DATA)),
                  SWE = double(NROW(DATA)),
                  NAOK=FALSE, DUP=FALSE, PACKAGE="hydromad")
        Sdischarge <- ans$U
        SWE <- ans$SWE
        ## make it a time series object again
        mostattributes(Sdischarge) <- attributes(DATA)
        class(Sdischarge) <- "ts"
        attributes(SWE) <- attributes(Sdischarge)
    } else {
        ## implementation in R for cross-checking (slow)
        ## time loop
        SWE <- Sdischarge <- P * 0
        LSWEprev <- LSWE_0
        ISWEprev <- ISWE_0
        for (t in seq(1, length(P))) {
            ## Melt (degree day model)
            melt <- min(max(kd*(E[t]-Tmelt),0),ISWEprev)
            ##Freezing (degree day model)
            freeze <- min(max(kf*(Tmelt-E[t]),0),LSWEprev)
            ##Mass balance for the snowpack
            ##
            ## Ice in the snowpack
            ISWE <- ISWEprev+Psnow[t]+freeze-melt
            ## Water in the snowpack
            LSWE <- min(rcap*ISWE,LSWEprev+Prain[t]+melt-freeze)
            ##Rain/melt is snowmelt discharge when there is snow on the ground,
            ##and rainfall in snow-free periods.
            Sdischarge[t] <- max(Prain[t]+melt-freeze-(rcap*ISWE-LSWEprev),0)
            SWE[t] <- LSWE+ISWE
            ISWEprev <- ISWE
            LSWEprev <- LSWE
        }
    }

    DATA[,"P"] <- Sdischarge

    ## IHACRES CMD-module
    if (return_state) {
        U <- cmd.sim(DATA, d = d, f = f, e = e, M_0 = M_0, return_state = TRUE)
        ans <- ts.union(U, SWE = SWE, TF = Sdischarge)
        colnames(ans)[1:NCOL(U)] <- colnames(U)
        return(ans)
    } else {
        U <- cmd.sim(DATA, d = d, f = f, e = e, M_0 = M_0, return_state = FALSE)
        return(U)
    }
}

snow.ranges <- function()
    list(Tmax=c(0, 2),
         Tmin=c(-1, 1),
         cr=c(0.8, 2),
         cs=c(0.8, 2),
         kd=c(2, 5),
         kf=c(0, 2),
         rcap=c(0, 1),
         d = c(50, 550),
         f = c(0.01, 3),
         e = c(0.01, 1.5))
