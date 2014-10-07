##Expects columns:
## Catchment calib.period sim.period Model.str     + objectives + other
## Objectives need to be specified (default r.squared)
## Assume objectives higher value is better. Values should be transformed prior to use
areModelsDominated <- function(res,objectives="r.squared"){
    res <- as.data.frame(res)
    
    if(!all(objectives %in% names(res))) stop(sprintf("objectives not recognised: %s",
                                                      paste(setdiff(objectives,names(res))),collapse=", "))
      
    id.vars <- intersect(names(res),c("Model.str","Catchment","calib.period","sim.period","Cal.objfn"))
    res <- res[,c(id.vars,objectives)]

    resm <- melt(res,
                 id.vars=id.vars,
                 measure.vars=objectives,
                 variable_name="objective"
                 )
    res2 <- cast(resm,...~sim.period)
    sim.periods <- unique(resm$sim.period)

    res2$dominated <- !paretoFilter(-as.matrix(res2[,sim.periods]))
    res2
}
