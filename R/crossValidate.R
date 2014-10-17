## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <josephguillaume@gmail.com>
##

crossValidate <- function(MODEL,periods,
                          name.Model.str=paste(MODEL$sma,MODEL$routing),
                          name.Cal.objfn="unknown",
                          name.Catchment=as.character(MODEL$call$DATA),
                          fitBy,...,
                          parallel=hydromad.getOption("parallel")[["crossValidate"]],
                          async=FALSE,export=c()
){
  cv.set <- splitData(MODEL,periods=periods)
  switch(parallel,
         "foreach"={
           opts=hydromad.options()
           runs=foreach(n=names(cv.set),
                        .packages="hydromad",
                        .inorder=FALSE,
                        .export=export,
                        .final=function(runs){
                          runs <- do.call(c,runs)
                          class(runs) <- unique(c("crossvalidation",class(runs),"runlist", "list"))
                          runs
                        },
                        .options.redis=list(async=async)
           ) %dopar% {
             hydromad.options(opts)
             if(isTRUE(hydromad.getOption("trace"))) cat("\nFitting period: ",n,"\n")    
             fitx <- fitBy(cv.set[[n]],...)
             new.runs <- update(cv.set,newpars=coef(fitx))
             names(new.runs) <- sprintf("%s_cal%s",names(cv.set),n)
             ## Preserve fit attributes
             new.runs[[sprintf("%s_cal%s",n,n)]] <- fitx
             for(m in 1:length(new.runs)){
               new.runs[[m]]$name.Model.str <- name.Model.str
               new.runs[[m]]$name.Cal.objfn <- name.Cal.objfn
               new.runs[[m]]$name.calib.period <- n
               new.runs[[m]]$name.sim.period <- names(cv.set)[m]
               new.runs[[m]]$name.Catchment <- name.Catchment
             }
             new.runs
           }
         },
{
  runs <- runlist()                   
  for(n in names(cv.set)){
    if(isTRUE(hydromad.getOption("trace"))) cat("\nFitting period: ",n,"\n")    
    fitx <- fitBy(cv.set[[n]],...)
    new.runs <- update(cv.set,newpars=coef(fitx))
    names(new.runs) <- sprintf("%s_cal%s",names(cv.set),n)
    ## Preserve fit attributes
    new.runs[[sprintf("%s_cal%s",n,n)]] <- fitx
    for(m in 1:length(new.runs)){
      new.runs[[m]]$name.Model.str <- name.Model.str
      new.runs[[m]]$name.Cal.objfn <- name.Cal.objfn
      new.runs[[m]]$name.calib.period <- n
      new.runs[[m]]$name.sim.period <- names(cv.set)[m]
      new.runs[[m]]$name.Catchment <- name.Catchment
    }
    runs <- c(runs,new.runs)
  }
  class(runs) <- unique(c("crossvalidation",class(runs),"runlist", "list"))
}
  ) ## switch parallel
return(runs)
}

summary.crossvalidation <- function(object, ...){
  s=NextMethod(object,...)
  s$sim.period <- sapply(object,function(x) x$name.sim.period)
  s$calib.period <- sapply(object,function(x) x$name.calib.period)
  s$Model.str <- sapply(object,function(x) x$name.Model.str)
  s$Cal.objfn <- sapply(object,function(x) x$name.Cal.objfn)
  s$Catchment <- sapply(object,function(x) x$name.Catchment)
  s
}  
