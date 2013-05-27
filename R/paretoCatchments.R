paretoCatchments <- function(stat,csv.name=NA,objectives="r.squared"){
    stat <- as.data.frame(stat)
    stat.split <- split(stat,stat$Catchment)
    stat.split <- lapply(stat.split,areModelsDominated,objectives=objectives)

    cat("\nHow many dominated\n")
    print(sapply(stat.split,function(res) table(res$dominated)))

    if (!is.na(csv.name)){
        cat("\nWrite csv of model realisations and whether dominated for each catchment\n")
        isd <- do.call(rbind,stat.split)
        isd$calib.period <- period.labels[as.character(isd$calib.period)]
        isd$dominated <- ifelse(isd$dominated,"Yes","No")
        isd <- isd[with(isd,order(Catchment,Model.str,calib.period)),]
        write.csv(isd,sprintf("%s_isdominated_models_catchments.csv",csv.name),row.names=FALSE)
    }

    cat("\nNon dominated model structures in each catchment\n")
    ##model.ok <- do.call(rbind,lapply(stat.split,function(xx) sort(unique(stat$Model.str)) %in% xx$Model.str[!xx$dominated]))
    model.ok <- do.call(rbind,lapply(stat.split,function(xx) {
        xx$Model.str <- factor(xx$Model.str)
        table(xx$Model.str[!xx$dominated])
    }))
    colnames(model.ok) <- sort(unique(stat$Model.str))
    ##model.ok <- ifelse(model.ok,"x","")
    print(model.ok)

    cat("\nNondominated models for each catchment\n")
    xx <- lapply(stat.split,function(res.rsq) res.rsq[!res.rsq$dominated,
                                                      setdiff(names(res.rsq),
                                                              c("objective","dominated"))
                                                      ])
    xx <- do.call(rbind,xx)
    if(!is.na(csv.name)) write.csv(xx,sprintf("%s_nondominated_models_catchments.csv",csv.name),row.names=FALSE)
    print(xx)

    cat("\nMinimum performance in each period\n")
    nondom.long <- do.call("rbind",lapply(stat.split,getDominatedLong))
    stopifnot(!any(is.na(nondom.long$dominated)))
    nondom.long <- nondom.long[!nondom.long$dominated,]
    min.perf.periods <- merge(
                              aggregate(value~sim.period+Catchment,data=nondom.long,min),
                              aggregate(value~sim.period+Catchment,data=nondom.long,max),
                              by=c("Catchment","sim.period"),suffixes=c(".min",".max")
                              )
    ## FIXME min.perf.periods$sim.period <- period.labels[as.character(min.perf.periods$sim.period)]
    min.perf.periods <- min.perf.periods[with(min.perf.periods,
                                              order(Catchment,sim.period)),]
    print(min.perf.periods)

    cat("\nMinimum performance across each period\n")
    min.perf.all <- merge(
                          aggregate(value~Catchment,data=nondom.long,min),
                          aggregate(value~Catchment,data=nondom.long,max),
                          by=c("Catchment"),suffixes=c(".min",".max")
                          )
    min.perf.all$range <- min.perf.all$value.max-min.perf.all$value.min
    print(min.perf.all)
    invisible(NULL)
}

getDominatedLong <- function(res.dom){
  resm2 <- as.data.frame(res.dom)
  resm2 <- melt(resm2,
                id.vars=intersect(names(res.dom),c("Model.str", "Catchment", "calib.period", "Cal.objfn","objective","dominated")),
                variable_name="sim.period")
  ## FIXME: ordering of sim.periods
  ##resm2$sim.period <- ordered(resm2$sim.period,levels=c("70","80","90","00"))
  ## FIXME: naming scheme of model realisations
  ##resm2$model <- paste(resm2$Model.str,resm2$calib.period,resm2$Cal.objfn)
  resm2
}
