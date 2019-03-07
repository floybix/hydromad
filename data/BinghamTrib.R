
BinghamTrib <- local({

    ## read files from doc directory
    pqdat <- utils::read.table("BinghamTrib.dat", header = TRUE,
                        col.names = c("Date", "P", "Q"), as.is = TRUE)
    tdat <- utils::read.table("DonnybrookTemp.txt", header = TRUE,
                       col.names = c("Year", "Month", "T"))
    ## convert dates
    pqdat$Date <- as.Date(pqdat$Date)
    tdat$Date <- with(tdat, as.Date(paste(Year, Month, 1, sep = "-")))

    ## convert from m^3/day to mm/day (i.e. divide by 1000 & catchment area km^2)
    pqdat$Q <- (pqdat$Q / 1000) / 2.68
#    pqdat$Q <- convertFlow(pqdat$Q, from = "m^3", area.km2 = 2.68)

    ## zoo objects
    library(zoo)
    tsPQ <- zoo::zoo(pqdat[,2:3], order.by = pqdat$Date, frequency = 1)
    tsT <- zoo::zoo(tdat$T, order.by = tdat$Date, frequency = 1)

    ## start temperature series from first month of PQ series
    tsT <- stats::window(tsT, start = zoo:::as.Date.yearmon(zoo::as.yearmon(stats::start(tsPQ))))
    ## make monthly temperature record a regular series (expand gaps)
    tsT <- merge(tsT, zoo::zoo(order.by = seq(stats::start(tsT), stats::end(tsT), by = "months")))
    ## fill monthly temperature gaps with seasonal averages
    months <- months(stats::time(tsT))
    avgs <- tapply(zoo::coredata(tsT), months, mean, na.rm = TRUE)
    tsT[is.na(tsT)] <- avgs[ match(months[is.na(tsT)], names(avgs)) ]

    tsPQE <- merge(tsPQ, E = tsT, all = TRUE)
    tsPQE$E <- zoo::na.locf(tsPQE$E, na.rm = FALSE)
    zoo::na.trim(tsPQE)
})
